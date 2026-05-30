# PharmCAT_annotation

Repository containing a script for merging, liftover, and annotation using the **PharmCAT** tool.

bash pharm3_new.sh /absolute/path/to/vcf_folder

This script merges individual `vcf.gz` files with a common population VCF file, performs normalization and liftover to **hg38** (script using Picard is pharm3_new.sh; script using CrossMap is pharm3_new_crossm.sh), identifies variants present in the **PharmCAT** database, annotates pharmacogene haplotypes, and determines phenotypes based on those haplotypes using **PharmCAT (Pharmacogenomics Clinical Annotation Tool)**. Generated files are saved in the `output` directory (created inside the folder containing the input `vcf.gz` files).

The pipeline performs the following steps:

- The script normalizes (saving statistics to `normalization_stats.tsv`), indexes, and merges individual `vcf.gz` files with the common population VCF (`vcf_jews_bgz.vcf.gz`);
- Positions listed in `stop_list.tsv` are removed from the merged VCF. The result is saved as `merged_excl_stop.vcf`;
- The filtered file is normalized again against the hg19 reference genome (`~/reference/hg19/hg19.fa`) and saved as `merged_excl_stop_normalized.vcf.gz`;
- CrossMap is used to convert variant coordinates from hg19 to hg38 using `hg19ToHg38.over.chain.gz`. Variants that cannot be converted are saved in `rejected_liftover38.vcf.gz`;
- PharmCAT preprocessing is performed using `pharmcat_vcf_preprocessor`, creating a separate VCF containing variants present in the PharmCAT database. Preprocessing is run in two modes:
  1. With the `-ss` (single sample) flag — separate VCF files are generated for each sample and saved in `output/pharmcat_preprocessed_vcfs`. File names correspond to sample IDs (used later for generating per-sample JSON reports);
  2. Without sample splitting — all PharmCAT variants are saved in a single VCF file in `output/pharmcat_preprocessed_all` (used for counting the total number of variants found in the PharmCAT database);
- PharmCAT JSON reports are generated. For each preprocessed VCF file, `pharmcat.jar` is executed to produce:
  1. Intermediate report `*.match.json` — contains the list of diplotypes identified in the VCF file for each pharmacogene;
  2. Phenotype report `*.phenotype.json` — assigns phenotypes based on diplotypes from `.match.json`;
  3. Final report `*.report.json` — maps phenotypes to clinical guideline databases (CPIC, DPWG) and contains diplotypes, phenotypes, and drug prescribing recommendations;
- JSON-to-TSV parsing. In the JSON report directory (`output/pharmcat_json_reports`), the Python script `pharm_script_new.py` extracts data from `*.report.json` files;
- Results are saved to `output/pharmcat_json_reports/tsv_output`:
  - `wide_diplotypes_phenotypes.tsv` — wide-format summary table (one row per sample) containing diplotypes, phenotypes, associated drugs, and SNPs for each gene;
    <img width="745" height="340" alt="568543778-0d173f40-0a75-439b-a1e7-ff2184310892" src="https://github.com/user-attachments/assets/1ed43a76-a951-4e1c-aba8-f0cfe46e20d2" />

  - `drugs_recommendations.tsv` — table with drug recommendations based on predicted phenotypes.
  - <img width="1401" height="479" alt="568550924-d2fcce7d-04be-49fe-bda5-3120d85aa482" src="https://github.com/user-attachments/assets/d9bd72e6-d1b5-47f1-9a0c-a407008b2664" />


---

All major pipeline steps, as well as variant counts at every processing stage, are recorded in `pipeline.log`.

Normalization statistics for each sample are stored in `normalization_stats.tsv`.

---

## Software Requirements

- Java 17 or newer
- bcftools 1.20
- htslib >= 1.18
- CrossMap
- PharmCAT 3.2.0
- Python 3.10.14 or newer

Required Python libraries:

- json
- glob
- csv
- os
- sys
- collections

---

## Reference Files

The script expects reference genomes and chain files in the following locations:

```bash
~/reference/hg19/hg19.fa
~/reference/hg38/hg38.fa
~/reference/hg38/hg19ToHg38.over.chain.gz
```

### Reference files must be indexed

- hg19.fa, hg19.fa.fai, hg19.dict
- hg38.fa, hg38.fa.fai, hg38.dict

---

## Picard Execution

Picard is executed within the script as follows:

```bash
java -jar ~/picard/build/libs/picard.jar LiftoverVcf \
    I="$NORMALIZED" \
    O="$LIFT_OVERED38" \
    CHAIN=~/reference/hg38/hg19ToHg38.over.chain.gz \
    REJECT="$REJECTED_LIFTOVER38" \
    R=~/reference/hg38/hg38.fa 2>> "$LOG_FILE"
```

---

## PharmCAT Installation

PharmCAT must be installed beforehand according to the official instructions:

```bash
curl -fsSL https://get.pharmcat.org | bash
```

```bash
pip3 install -r requirements.txt
```

---

## Running the Script

```bash
bash pharm3_new.sh ~/absolute/path
```

where `~/absolute/path` is the absolute path to the directory containing the VCF files to be processed.

---

## Required Files in the Working Directory

The directory from which the script is executed must contain:

- Input VCF files to be processed
- Population VCF file (`vcf_jews_bgz.vcf.gz`)
- Stop-list variants table (`stop_list.tsv`)
- Pipeline script (`pharm3_new.sh`)
- Python script for processing PharmCAT JSON reports
<img width="680" height="364" alt="568561002-024231e7-a8cf-4eac-931b-2329fda8feab" src="https://github.com/user-attachments/assets/f03e9c96-e572-4fc7-ae02-ee5f119c9433" />

---

## Output Files

After the script finishes, the `output` directory (created inside the sample directory) will contain:

### VCF Files

- `merged.vcf.gz` — merged VCF containing all samples (unfiltered)
- `merged_excl_stop.vcf.gz` — merged VCF after removing stop-list variants
- `merged_excl_stop_normalized.vcf.gz` — normalized VCF
- `merged_excl_stop.vcf.gz.csi` — index file
- `lift38_merged_excl_stop_normalized.vcf.gz` — liftover result (hg19 → hg38)
- `lift38_merged_excl_stop_normalized.vcf.gz.tbi` — index file
- `lift38_rejected.vcf.gz` — variants rejected during liftover
- <img width="704" height="503" alt="565552430-66fa645d-487f-41bd-bf2a-f65ed847cd8e" src="https://github.com/user-attachments/assets/874798b7-dec6-4ef1-9513-77e98d68bc3e" />


### Logs and Statistics

- `pipeline.log` — pipeline log file
- `normalization_stats.tsv` — normalization statistics table containing:
  - Sample
  - Total
  - Split
  - Joined
  - Realigned
  - Skipped

### PharmCAT Outputs

- `pharmcat_preprocessed_vcfs/` — sample-specific preprocessed VCF files containing variants found in the PharmCAT database
- `pharmcat_preprocessed_all/` — combined VCF file containing all samples and PharmCAT variants
- `pharmcat_json_reports/` — PharmCAT JSON reports

---

## Final Results

Final output tables are located in:

```text
output/pharmcat_json_reports/tsv_output
```

This directory contains:

- `wide_diplotypes_phenotypes.tsv` — summary table with diplotypes, phenotypes, associated drugs, and SNPs
- `drugs_recommendations.tsv` — phenotype-based drug recommendation table
- <img width="326" height="113" alt="568552927-b528f11a-8778-4026-a587-5ca2da15ed2f" src="https://github.com/user-attachments/assets/927d0a92-0d40-4455-82e4-d94d5e3be3c1" />


These files represent the final pharmacogenomic interpretation generated from the input VCF datasets.
