#!/bin/bash
set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <vcf_folder> [sample_to_remove]"
    exit 1
fi

VCF_DIR="$1"
OUT_DIR="$VCF_DIR/output"
PROC_DIR="$OUT_DIR/pharmcat_preprocessed_vcfs"
PROC_ALL_DIR="$OUT_DIR/pharmcat_preprocessed_all"
JSON_DIR="$OUT_DIR/pharmcat_json_reports"

mkdir -p "$OUT_DIR" "$PROC_DIR" "$PROC_ALL_DIR" "$JSON_DIR"

if [ -f "$VCF_DIR/pharm_script_new.py" ]; then
    cp "$VCF_DIR/pharm_script_new.py" "$JSON_DIR/"
else
    echo "WARNING: pharm_script_new.py not found in $VCF_DIR"
fi

MERGED="$OUT_DIR/merged.vcf.gz"
MERGED_EXCL_STOP="$OUT_DIR/merged_excl_stop.vcf.gz"
NORMALIZED="$OUT_DIR/merged_excl_stop_normalized.vcf.gz"
LIFT_OVERED38="$OUT_DIR/lift38_merged_excl_stop_normalized.vcf.gz"
REJECTED_LIFTOVER38="$OUT_DIR/lift38_rejected.vcf.gz"
SAMPLES="$PROC_DIR/input_samples.txt"
SCRIPT_JSON="$JSON_DIR/pharm_script_new.py"

STATS_FILE="$OUT_DIR/normalization_stats.tsv"
echo -e "Sample\tTotal\tSplit\tJoined\tRealigned\tSkipped" > "$STATS_FILE"

LOG_FILE="$OUT_DIR/pipeline.log"

echo "STEP 1: SORT, INDEX AND NORMALIZE ALL VCFs"

shopt -s nullglob
files=( "$VCF_DIR"/*.vcf.gz )

if [ ${#files[@]} -eq 0 ]; then
    echo "ERROR: No .vcf.gz files found in $VCF_DIR"
    exit 1
fi

for f in "${files[@]}"; do
    base=$(basename "$f" .vcf.gz)
    echo "→ Processing $base"

    norm_file="$OUT_DIR/${base}.sorted_norm.vcf.gz"
    sorted_file="$OUT_DIR/${base}.sorted.vcf.gz"

    bcftools sort "$f" -Oz -o "$sorted_file" 2>> "$LOG_FILE"
    norm_output=$(bcftools norm "$sorted_file" -m-any -Oz -o "$norm_file" 2>&1)

    rm -f "$sorted_file"
    bcftools index "$norm_file"

    line=$(echo "$norm_output" | grep "Lines.*total/split/joined/realigned/skipped" | tail -n1)

    numbers=$(echo "$line" | sed -E 's/.*total\/split\/joined\/realigned\/skipped:[[:space:]]*([0-9\/]+).*/\1/')

    if [ -z "$numbers" ]; then
        total=0; split=0; joined=0; realigned=0; skipped=0
    else
        IFS='/' read -r total split joined realigned skipped <<< "$numbers"
    fi

    echo -e "${base}\t${total}\t${split}\t${joined}\t${realigned}\t${skipped}" >> "$STATS_FILE"
done

echo "STEP 2: MERGE of VCFs"
bcftools merge "$OUT_DIR"/*.sorted_norm.vcf.gz -Oz -o "$MERGED"
bcftools index "$MERGED"

echo "Variants number after merge of VCFs: $(bcftools view -H "$MERGED" | wc -l)" >> "$LOG_FILE"

rm -f "$OUT_DIR"/*.sorted_norm.vcf.gz "$OUT_DIR"/*.csi

echo "STEP 3: FILTER OUT STOP POSITIONS"

bcftools view -H "$MERGED" | \
awk 'NR==FNR{a[$1":"$2":"$3":"$4]; next} !($1":"$2":"$4":"$5 in a)' stop_list.tsv - \
| cat <(bcftools view -h "$MERGED") - | bgzip -c > "$MERGED_EXCL_STOP"

echo "Variants number after filtering out of stop positions: $(bcftools view -H  "$MERGED_EXCL_STOP" | wc -l)" >> "$LOG_FILE"

echo "STEP 4: NORMALIZATION"

bcftools norm -m-any -f ~/reference/hg19/hg19.fa \
    -Oz -o "$NORMALIZED" "$MERGED_EXCL_STOP" 2>> "$LOG_FILE"

bcftools index "$NORMALIZED"

echo "Variants number after normalization: $(bcftools view -H "$NORMALIZED" | wc -l)" >> "$LOG_FILE"

echo "STEP 5: LIFT OVER TO HG38"

java -jar ~/picard/build/libs/picard.jar LiftoverVcf \
    I="$NORMALIZED" \
    O="$LIFT_OVERED38" \
    CHAIN=~/reference/hg38/hg19ToHg38.over.chain.gz \
    REJECT="$REJECTED_LIFTOVER38" \
    R=~/reference/hg38/hg38.fa 2>> "$LOG_FILE"

echo "Variants number after liftover to hg38: $(bcftools view -H "$LIFT_OVERED38" | wc -l)" >> "$LOG_FILE"

echo "STEP 6: PHARMCAT PREPROCESSING"

~/pharmcat/pharmcat_vcf_preprocessor \
    -vcf "$LIFT_OVERED38" \
    -ss \
    --output-dir "$PROC_DIR" 2>> "$LOG_FILE"

ls "$PROC_DIR"/*.preprocessed.vcf | xargs -n1 basename > "$SAMPLES"

~/pharmcat/pharmcat_vcf_preprocessor -vcf "$LIFT_OVERED38" --output-dir "$PROC_ALL_DIR" 2>> "$LOG_FILE"

echo "Variants in PharmCAT preprocessed VCF: $(grep -h -v "^#" "$PROC_ALL_DIR"/*.vcf | wc -l)" >> "$LOG_FILE"

echo "STEP 7: PHARMCAT JSON REPORTS"

cd "$PROC_DIR"
while read -r i; do
    SAMPLE=${i%.preprocessed.vcf}

    java -jar ~/pharmcat/pharmcat.jar \
        -vcf "$PROC_DIR/$i" \
        -bf "$SAMPLE" \
        -reporterJson \
        -o "$JSON_DIR" 2>> "$LOG_FILE"

done < "$SAMPLES"

echo "STEP 8: PARSE JSON REPORTS"

cd "$JSON_DIR"
python "$SCRIPT_JSON" >> "$JSON_DIR/parse_json_log.txt" 2>&1

echo "DONE → Results (tsv tables) are in $JSON_DIR/tsv_output"