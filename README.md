# PharmCAT_annotation
Repository with script for merge, liftover and annotation with PharmCAT tool

Скрипт для объединения vcf.gz файлов с общим vcf.gz файлом с образцами, нормализации и liftover  до hg38, поиска вариантов из базы PharmCAT, аннотации гаплотипов фармакогенов (и опредления фенотипов согласно гаплотипам) при помощи PharmCAT (Pharmacogenomics Clinical Annotation Tool). Полученные файлы сохраняются в папку output (создается внутри папки с vcf.gz файлами) <br>
<br>
- скрипт нормализует (данные по нормализации сохраняет в normalization_stats.tsv), индексирует и объединяет отдельные vcf.gz файлы с общим vcf.gz файлом (vcf_jews_bgz.vcf.gz);  <br>
- из объединённого vcf удаляются позиции, перечисленные в файле stop_list.tsv, результат сохраняется как merged_excl_stop.vcf; <br>
- отфильтрованный файл снова нормализуется с привязкой к референсному геному hg19 (~/reference/hg19/hg19.fa) и сохраняется, как merged_excl_stop_normalized.vcf.gz; <br>
- с помощью CrossMap выполняется конвертация координат вариантов из hg19 в hg38, используется hg19ToHg38.over.chain.gz. Варианты, которые не удалось конвертировать, сохраняются в отдельный файл rejected_liftover38.vcf.gz; <br>
- предобработка для PharmCAT с помощью pharmcat_vcf_preprocessor - в отдельный vcf  файл сохраняются варианты, присутствующие в базе данных PharmCAT. Предобработка проводится в двух режимах: <br>
  1) с флагом -ss (single sample) — результат сохраняется в виде отдельных vcf файлов по образцам в папку output/pharmcat_preprocessed_vcfs, имена файлов соответствуют образцам (это делается для будущей генерации JSON отчетов по образцам); <br>
  2) также проводится предобработка без разделения по образцам — все варианты сохраняются в один общий VCF  файл в папку output/pharmcat_preprocessed_all (это делается для подсчета общего числа вариантов, обнаруженных в базе PharmCAT); <br>
  <br>
- генерация JSON-отчётов PharmCAT. Для каждого предобработанного VCF-файла запускается pharmcat.jar для генерации отчетов (сохраняются в папку output/pharmcat_json_reports): <br>
  1) промежуточный отчет *.match.json - содержит список диплотипов, которые были найдены в VCF-файле для каждого фармакогена; <br>
  2) отчет с результатами фенотпирования (прогнозирование фенотипа) *.phenotype.json  - берет диплотипы из .match.json и присваивает им  фенотипы; <br>
  3) итоговый отчет *.report.json - берет фенотипы из .phenotype.json и сопоставляет их с базами данных клинических рекомендаций (CPIC, DPWG), содержит не только диплотипы и фенотипы, но и рекомендации по выбору препарата. <br>
- парсинг JSON в TSV. В папке с JSON-отчётами (output/pharmcat_json_reports) запускается Python-скрипт pharm_script_new.py, которые извлекает данные из файлов *.report.json
- результаты сохраняются в подпапке output/pharmcat_json_reports/tsv_output: wide_diplotypes_phenotypes.tsv — сводная таблица в "широком" формате (одна строка на образец, колонки для каждого гена с диплотипом, фенотипом, лекарственными препаратами и SNP), drugs_recommendations.tsv - таблица с рекомендациями по приему лекарственных препаратов в зависимости от фенотипа

<img width="745" height="340" alt="image" src="https://github.com/user-attachments/assets/0d173f40-0a75-439b-a1e7-ff2184310892" />

<br>

<img width="1401" height="479" alt="image" src="https://github.com/user-attachments/assets/d2fcce7d-04be-49fe-bda5-3120d85aa482" />


<br>
<br>
Все ключевые шаги, а также подсчет количества вариантов в vcf файлах на кждом этапе записываются в pipeline.log. Для каждого образца сохраняется статистика нормализации в normalization_stats.tsv
<br>
<br>

**Программы, которые используются при выполнении скрипта: <br>**
- Java 17 или новее
- bcftools 1.20
- htslib >= v1.18
- CrossMap
- PharmCAT 3.2.0
- Python  3.10.14 или новее (библиотеки json, glob, csv, os, sys, collections) <br>

**В скрипте доступ к референсным последовательностям и chain, необходимыем для нормализации и liftover, проводится таким образом:**
```bash
~/reference/hg19/hg19.fa
~/reference/hg38/hg38.fa
~/reference/hg38/hg19ToHg38.over.chain.gz
```
**Референсные файлы должны быть индексированы:**
- hg19.fa, hg19.fa.fai, hg19.dict
- hg38.fa, hg38.fa.fai, hg38.dict

**Picard в скрипте запускается таким образом:**
```bash
java -jar ~/picard/build/libs/picard.jar LiftoverVcf \
    I="$NORMALIZED" \
    O="$LIFT_OVERED38" \
    CHAIN=~/reference/hg38/hg19ToHg38.over.chain.gz \
    REJECT="$REJECTED_LIFTOVER38" \
    R=~/reference/hg38/hg38.fa 2>> "$LOG_FILE"
```

**PharmCAT необходимо предварительно установить при помощи (https://pharmcat.clinpgx.org/using/Setup-PharmCAT/): <br>**
```bash
curl -fsSL https://get.pharmcat.org | bash
```
```python
pip3 install -r requirements.txt
```

**Запуск скрипта:** bash pharm3_new.sh  ~/absolute/path (абсолютный путь к папке с vcf файлами, которые нужно обработать) <br>
 <br>
**В папке, в которой запускается скрипт,  находятся:** <br>
- vcf файлы, которые нужно обработать <br>
- общий vcf файл по популяциям (vcf_jews_bgz.vcf.gz)  <br>
- таблица с вариантами из стоп-листа, которые нужно исключить (stop_list.tsv) <br>
- скрипт pharm3_new.sh   <br>
- скрипт для обработки JSON-отчетов (.py)
<img width="680" height="364" alt="image" src="https://github.com/user-attachments/assets/024231e7-a8cf-4eac-931b-2329fda8feab" />





**После выполнения скрипта в папке output (будет создана в папке с образцами) будут находиться:** <br>
- объединенный файл со всеми образцами с индексацией (файл без фильтрации, исходный после объединения) - merged.vcf.gz
- объединенный файл после исключения вариантов из стоп-листа - merged_excl_stop.vcf.gz
- нормализованный merged_excl_stop.vcf.gz - merged_excl_stop_normalized.vcf.gz и его индекс merged_excl_stop.vcf.gz.csi
- файл после liftover с hg19 до hg 38: lift38_merged_excl_stop_normalized.vcf.gz и его индекс lift38_merged_excl_stop_normalized.vcf.gz.tbi
- файл с вариантами, которые не удалось конвертировать (были отклонены при liftover) - lift38_rejected.vcf.gz
- лог файл pipeline.log
- таблица с данными при нормализации vcf файлов (Sample, Total, Split, Joined, Realigned,	Skipped) -  normalization_stats.tsv
- папка pharmcat_preprocessed_vcfs с отдельными подготовленными vcf файлами по образцам (получены при помощи pharmcat_vcf_preprocessor, содержат варианты, обнаруженные в базе PharmCAT)
- папка pharmcat_preprocessed_all с общим vcf файлом со всеми образцами (содержит варианты, обнаруженные в базе PharmCAT)
- папка с JSON-отчетами PharmCAT - pharmcat_json_reports

  <img width="704" height="503" alt="image" src="https://github.com/user-attachments/assets/66fa645d-487f-41bd-bf2a-f65ed847cd8e" />


**Результат (таблицы all_genes_combined.tsv и samples_by_gene_pivot.tsv) содержится в папке output/pharmcat_json_reports/tsv_output**
<img width="326" height="113" alt="image" src="https://github.com/user-attachments/assets/b528f11a-8778-4026-a587-5ca2da15ed2f" />

