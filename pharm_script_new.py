#!/usr/bin/env python3
import json
import pandas as pd
import re
import os
import sys  # Добавлен импорт sys
from pathlib import Path

# Исправлена логика путей
input_folder = sys.argv[1] if len(sys.argv) > 1 else "."
output_dir = os.path.join(input_folder, "tsv_output") # Исправлено на input_folder
os.makedirs(output_dir, exist_ok=True)

all_data_list = []
all_drugs_list = []

if not os.path.exists(input_folder):
    print(f"Папка {input_folder} не найдена!")
    sys.exit(1)

# Собираем файлы
files = list(Path(input_folder).glob('*.report.json'))
if not files:
    print("Файлы *.report.json не найдены.")
    sys.exit(0)

for file_path in files:
    with open(file_path, 'r', encoding='utf-8') as f:
        try:
            data = json.load(f)
        except Exception as e:
            print(f"Ошибка в файле {file_path}: {e}")
            continue
    
    match = re.search(r'\.([a-z0-9]+)\.', file_path.name)
    s_id = match.group(1) if match else file_path.stem

    sample_phenotypes = {} 
    genes_dict = data.get('genes', {})
    
    for g_symbol, g_info in genes_dict.items():
        sd_list = g_info.get('sourceDiplotypes', [])
        diplo = sd_list[0] if isinstance(sd_list, list) and len(sd_list) > 0 else {}
        
        p_list = diplo.get('phenotypes', [])
        p_val = p_list[0] if p_list else "Unknown"
        sample_phenotypes[g_symbol] = p_val

        variants = g_info.get('variants', [])
        snps = [f"[{v.get('dbSnpId')}: {v.get('chromosome')}:{v.get('position')} {v.get('call')}]" for v in variants]
        
        all_data_list.append({
            'sample_id': s_id,
            'gene': g_symbol,
            'Diplotype': diplo.get('label', 'N/A'),
            'Phenotype': p_val,
            'RelatedDrugs': ", ".join([d.get('name', '') for d in g_info.get('relatedDrugs', [])]),
            'SNP_Details': "; ".join(snps) if snps else "No variants"
        })

    cpic_root = data.get('drugs', {}).get('CPIC Guideline Annotation', {})
    for drug_key, d_info in cpic_root.items():
        for guide in d_info.get('guidelines', []):
            for annot in guide.get('annotations', []):
                target_phenos = annot.get('phenotypes', {}) 
                if not target_phenos: continue
                
                is_match = True
                matched_parts = []
                for g_name, required_p in target_phenos.items():
                    current_p = sample_phenotypes.get(g_name)
                    matched_parts.append(f"{g_name}: {current_p}")
                    if current_p != required_p:
                        is_match = False
                        break
                
                if is_match:
                    all_drugs_list.append({
                        'sample_id': s_id,
                        'drug': d_info.get('name', drug_key),
                        'matched_phenotype': "; ".join(matched_parts),
                        'implications': " ".join(annot.get('implications', [])),
                        'recommendation': annot.get('drugRecommendation'),
                        'classification': annot.get('classification'),
                        'genes_involved': " + ".join(target_phenos.keys())
                    })

# Сохранение с проверкой на наличие данных
if all_data_list:
    df_raw = pd.DataFrame(all_data_list)
    df_pivot = df_raw.pivot_table(index='sample_id', columns='gene', 
                                  values=['Diplotype', 'Phenotype', 'RelatedDrugs', 'SNP_Details'], 
                                  aggfunc='first')
    df_pivot.columns = [f"{c[1]}_{c[0]}" for c in df_pivot.columns]
    df_pivot = df_pivot[sorted(df_pivot.columns)].reset_index()
    
    output_pivot = os.path.join(output_dir, 'wide_diplotypes_phenotypes.tsv')
    df_pivot.to_csv(output_pivot, sep='\t', index=False)
    print(f"Таблица фенотипов сохранена в: {output_pivot}")

if all_drugs_list:
    df_drugs = pd.DataFrame(all_drugs_list).drop_duplicates()
    output_drugs = os.path.join(output_dir, 'drugs_recommendations.tsv')
    cols_order = ['sample_id', 'drug', 'matched_phenotype', 'recommendation', 'implications', 'classification', 'genes_involved']
    df_drugs[cols_order].to_csv(output_drugs, sep='\t', index=False)
    print(f"Рекомендации по лекарствам сохранены в: {output_drugs}")
