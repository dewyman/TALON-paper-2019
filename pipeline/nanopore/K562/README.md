## Make a whitelist file of filtered GM12878 transcripts
```
python /pub/dwyman/TALON/post-TALON_tools/filter_talon_transcripts.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -p pairings.csv \
          --o K562_whitelist.csv
```

## Make a GTF file using the GM12878 whitelist
```
python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist K562_whitelist.csv \
          --datasets K562_datasets.txt \
          --o K562_filtered
```

## Make an abundance matrix without filtering
```
python /pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db ../full_gencode_v29_2019-05-24.db \
           -a gencode_v29 \
           --build hg38 \
           --o ../all
```

## Map antisense TALON genes to their sense counterparts
```
python /pub/dwyman/TALON/post-TALON_tools/map_antisense_genes_to_sense.py \
    --db ../full_gencode_v29_2019-05-24.db \
    --a gencode_v29 \
    --o ../all
```

## RNA-PET analysis
```
mkdir -p RNA-PET
module load bedtools
cd ../../../RNA-PET
python run_RNA-PET_analysis.py \
    --gtf ../pipeline/nanopore/K562/K562_filtered_talon.gtf \
    --rnapet data/K562_hg38.bed \
    --maxdist 100 \
    --o ../pipeline/nanopore/K562/RNA-PET/K562

cd ../pipeline/nanopore/K562
Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R  \
     --f RNA-PET/K562_RNA-PET_results.csv   \
     --t RNA-PET \
     --novelty RNA-PET/transcript_beds/K562_novelty.csv \
     --abundance ../all_talon_abundance.tsv \
     --d1 ONT31 --d2 ONT18 \
     --as ../all_antisense_mapping.csv \
     --splitISM \
     -o RNA-PET/K562
```

## CAGE analysis, ENCODE
```
source activate mypython3.7.2
mkdir -p CAGE/ENCODE/
python ../../../CAGE/run_CAGE_analysis.py \
        --gtf K562_filtered_talon.gtf \
        --cage ../../../CAGE/data/ENCODE/K562_CAGE.bed \
        --maxdist 100 \
        --o CAGE/ENCODE/K562

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R  \
     --f CAGE/ENCODE/K562_CAGE_results.csv   \
     --t CAGE \
     --novelty CAGE/ENCODE/transcript_beds/K562_novelty.csv \
     --abundance ../all_talon_abundance.tsv \
     --d1 ONT31 --d2 ONT18 \
     --as ../all_antisense_mapping.csv \
     --splitISM \
     -o CAGE/ENCODE/K562
```

## CAGE analysis, FANTOM
```
source activate mypython3.7.2
mkdir -p CAGE/FANTOM5/
python ../../../CAGE/run_CAGE_analysis.py \
        --gtf K562_filtered_talon.gtf \
        --cage ../../../CAGE/data/FANTOM5/hg38_CAGE.bed \
        --maxdist 100 \
        --o CAGE/FANTOM5/K562

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R  \
     --f CAGE/FANTOM5/K562_CAGE_results.csv   \
     --t CAGE \
     --novelty CAGE/FANTOM5/transcript_beds/K562_novelty.csv \
     --abundance ../all_talon_abundance.tsv \
     --d1 ONT31 --d2 ONT18 \
     --as ../all_antisense_mapping.csv \
     --splitISM \
     -o CAGE/FANTOM5/K562
```

## Computational PAS analysis
```
source activate mypython3.7.2
mkdir PAS-comp
python ../../../PAS-computational/run_computational_PAS_analysis.py \
        --gtf K562_filtered_talon.gtf \
        --genome ../../../refs/hg38/hg38.fa \
        --maxdist 35 \
        --o PAS-comp/K562

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R \
    --f PAS-comp/K562_polyA_motif.csv \
    --t PAS-comp \
    --novelty PAS-comp/transcript_beds/K562_novelty.csv \
    --abundance ../all_talon_abundance.tsv \
    --d1 ONT31 --d2 ONT18 --as ../all_antisense_mapping.csv \
    --splitISM \
    -o PAS-comp/K562
```

