#!/usr/bin/env python3
import json
import glob
import csv
import os
import sys
from collections import defaultdict

def load_json(json_file):
    with open(json_file, "r") as f:
        return json.load(f)

def get_sample_id_from_filename(filename):
    """
    Извлекает zy1874 из 'stop_filt_lifted_over38.zy1874.preprocessed.report.json'
    Если формат другой, возвращает имя файла без расширений.
    """
    basename = os.path.basename(filename)
    parts = basename.split('.')
    # Обычно ID идет вторым элементом в таком длинном имени
    if len(parts) > 2:
        return parts[1] 
    return basename.split('.')[0]

def extract_diplotype_info(gene_info):
    """Безопасное извлечение диплотипа, фенотипа и скора"""
    diplo_list = gene_info.get("sourceDiplotypes", [])
    diplo = diplo_list[0] if diplo_list else {}

    diplotype = diplo.get("label", "")
    
    # Безопасно берем первый фенотип из списка
    phenotypes = diplo.get("phenotypes", [])
    phenotype = phenotypes[0] if phenotypes else ""

    function = diplo.get("allele1", {}).get("function", "")
    
    # Конвертируем score в строку, избегая None
    activity = diplo.get("activityScore")
    activity = str(activity) if activity is not None else ""

    return diplotype, phenotype, function, activity

def extract_from_report(json_file):
    data = load_json(json_file)
    # Приоритет: ID из метаданных -> ID из имени файла
    sample_id = data.get("metadata", {}).get("sampleId") or get_sample_id_from_filename(json_file)
    
    rows = []
    for gene_name, gene_info in data.get("genes", {}).items():
        diplotype, phenotype, function, activity = extract_diplotype_info(gene_info)
        rows.append({
            "sample": sample_id,
            "gene": gene_name,
            "diplotype": diplotype,
            "phenotype": phenotype,
            "function": function,
            "activity_score": activity,
            "source_file": os.path.basename(json_file)
        })
    return rows

def extract_from_phenotype(json_file):
    data = load_json(json_file)
    sample_id = data.get("matcherMetadata", {}).get("sampleId") or get_sample_id_from_filename(json_file)
    
    rows = []
    for gene_name, gene_info in data.get("geneReports", {}).items():
        diplotype, phenotype, function, activity = extract_diplotype_info(gene_info)
        rows.append({
            "sample": sample_id, "gene": gene_name, "diplotype": diplotype,
            "phenotype": phenotype, "function": function, "activity_score": activity,
            "source_file": os.path.basename(json_file)
        })
    return rows

def extract_from_match(json_file):
    data = load_json(json_file)
    sample_id = data.get("metadata", {}).get("sampleId") or get_sample_id_from_filename(json_file)
    
    rows = []
    for result in data.get("results", []):
        gene = result.get("gene")
        if not gene: continue
        diplotypes = result.get("diplotypes", [])
        diplotype = diplotypes[0].get("name", "") if diplotypes else ""
        rows.append({
            "sample": sample_id, "gene": gene, "diplotype": diplotype,
            "phenotype": "", "function": "", "activity_score": "",
            "source_file": os.path.basename(json_file)
        })
    return rows

def main():
    input_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    output_dir = os.path.join(input_dir, "tsv_output")
    os.makedirs(output_dir, exist_ok=True)

    all_rows = []
    patterns = [
        ("*.report.json", extract_from_report),
        ("*.phenotype.json", extract_from_phenotype),
        ("*.match.json", extract_from_match),
    ]

    for pattern, extractor in patterns:
        files = glob.glob(os.path.join(input_dir, pattern))
        print(f"Found {len(files)} {pattern}")
        for file in files:
            try:
                rows = extractor(file)
                all_rows.extend(rows)
                print(f"  Processed {os.path.basename(file)}")
            except Exception as e:
                print(f"Error processing {file}: {e}", file=sys.stderr)

    if not all_rows:
        print("No data extracted!"); return

    # Дедупликация с учетом длинных имен файлов
    priority_map = {"report.json": 3, "phenotype.json": 2, "match.json": 1}
    unique = {}

    for row in all_rows:
        key = (row["sample"], row["gene"])
        fname = row["source_file"]
        
        # Определяем тип файла по расширению (надежный способ)
        this_type = "match.json"
        if "report.json" in fname: this_type = "report.json"
        elif "phenotype.json" in fname: this_type = "phenotype.json"

        old_row = unique.get(key)
        if not old_row:
            unique[key] = row
        else:
            old_fname = old_row["source_file"]
            old_type = "match.json"
            if "report.json" in old_fname: old_type = "report.json"
            elif "phenotype.json" in old_fname: old_type = "phenotype.json"
            
            if priority_map[this_type] > priority_map[old_type]:
                unique[key] = row

    unique_rows = sorted(unique.values(), key=lambda x: (x["sample"], x["gene"]))

    # Запись основного TSV
    combined_file = os.path.join(output_dir, "all_genes_combined.tsv")
    with open(combined_file, "w", newline="") as f:
        fieldnames = ["sample", "gene", "diplotype", "phenotype", "function", "activity_score", "source_file"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(unique_rows)

    # Пивот (одна строка на образец)
    samples = defaultdict(dict)
    for r in unique_rows:
        samples[r["sample"]][r["gene"]] = (r["diplotype"], r["phenotype"])

    genes = sorted({r["gene"] for r in unique_rows})
    pivot_file = os.path.join(output_dir, "samples_by_gene_pivot.tsv")
    with open(pivot_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        header = ["sample"]
        for g in genes: header += [f"{g}_diplotype", f"{g}_phenotype"]
        writer.writerow(header)
        for s in sorted(samples):
            row = [s]
            for g in genes:
                d, p = samples[s].get(g, ("Unknown", "Unknown"))
                row += [d, p]
            writer.writerow(row)

    print(f"\nSuccess! Files saved in: {output_dir}")

if __name__ == "__main__":
    main()
