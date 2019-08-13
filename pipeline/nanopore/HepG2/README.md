## Make a whitelist file of filtered HepG2 transcripts
```
python /pub/dwyman/TALON/post-TALON_tools/filter_talon_transcripts.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -p pairings.csv \
          --o HepG2_whitelist.csv
```

## Make a GTF file using the HepG2 whitelist
```
python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist HepG2_whitelist.csv \
          --datasets HepG2_datasets.txt \
          --o HepG2_filtered
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
    --gtf ../pipeline/nanopore/HepG2/HepG2_filtered_talon.gtf \
    --rnapet data/HepG2_hg38.bed \
    --maxdist 100 \
    --o ../pipeline/nanopore/HepG2/RNA-PET/HepG2

cd ../pipeline/nanopore/HepG2
Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R  \
     --f RNA-PET/HepG2_RNA-PET_results.csv   \
     --t RNA-PET \
     --novelty RNA-PET/transcript_beds/HepG2_novelty.csv \
     --abundance ../all_talon_abundance.tsv \
     --d1 ONT33 --d2 ONT34 \
     --as ../all_antisense_mapping.csv \
     --splitISM \
     -o RNA-PET/HepG2
```

## CAGE analysis, ENCODE
```
source activate mypython3.7.2
mkdir -p CAGE/ENCODE/
python ../../../CAGE/run_CAGE_analysis.py \
        --gtf HepG2_filtered_talon.gtf \
        --cage ../../../CAGE/data/ENCODE/HepG2_CAGE.bed \
        --maxdist 100 \
        --o CAGE/ENCODE/HepG2

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R  \
     --f CAGE/ENCODE/HepG2_CAGE_results.csv   \
     --t CAGE \
     --novelty CAGE/ENCODE/transcript_beds/HepG2_novelty.csv \
     --abundance ../all_talon_abundance.tsv \
     --d1 ONT33 --d2 ONT34 \
     --as ../all_antisense_mapping.csv \
     --splitISM \
     -o CAGE/ENCODE/HepG2
```

## CAGE analysis, FANTOM
```
source activate mypython3.7.2
mkdir -p CAGE/FANTOM5/
python ../../../CAGE/run_CAGE_analysis.py \
        --gtf HepG2_filtered_talon.gtf \
        --cage ../../../CAGE/data/FANTOM5/hg38_CAGE.bed \
        --maxdist 100 \
        --o CAGE/FANTOM5/HepG2

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R  \
     --f CAGE/FANTOM5/HepG2_CAGE_results.csv   \
     --t CAGE \
     --novelty CAGE/FANTOM5/transcript_beds/HepG2_novelty.csv \
     --abundance ../all_talon_abundance.tsv \
     --d1 ONT33 --d2 ONT34 \
     --as ../all_antisense_mapping.csv \
     --splitISM \
     -o CAGE/FANTOM5/HepG2
```

## Computational PAS analysis
```
source activate mypython3.7.2
mkdir PAS-comp
python ../../../PAS-computational/run_computational_PAS_analysis.py \
        --gtf HepG2_filtered_talon.gtf \
        --genome ../../../refs/hg38/hg38.fa \
        --maxdist 35 \
        --o PAS-comp/HepG2

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R \
    --f PAS-comp/HepG2_polyA_motif.csv \
    --t PAS-comp \
    --novelty PAS-comp/transcript_beds/HepG2_novelty.csv \
    --abundance ../all_talon_abundance.tsv \
    --d1 ONT33 --d2 ONT34 --as ../all_antisense_mapping.csv \
    --splitISM \
    -o PAS-comp/HepG2
```

## ONT vs. Illumina quantification analysis
```
# gene detection
Rscript ../../../compare_to_FLAIR/plot_detection_by_TPM_for_datasets.R \
    --f all_talon_abundance.tsv \
    --datasets ONT33,ONT32 \
    --ik1 ../../../Illumina/HepG2/Kallisto/Rep1/abundance.tsv \
    --ik2 ../../../Illumina/HepG2/Kallisto/Rep2/abundance.tsv \
    --color green \
    --outdir ../HepG2/

# differences in gene quantification
Rscript ../../../analysis_scripts/pacbio_v_illumina_edgeR.R \
    --f all_talon_abundance.tsv \
    --datasets ONT33,ONT32 \
    --ik1 ../../../Illumina/HepG2/Kallisto/Rep1/abundance.tsv \
    --ik2 ../../../Illumina/HepG2/Kallisto/Rep2/abundance.tsv \
    --dtype ONT \
    --color green \
    --outdir ../HepG2/

# differences in transcript quantification 

Rscript ../../../analysis_scripts/pacbio_v_illumina_edgeR.R \
    --f all_talon_abundance.tsv \
    --datasets ONT32,ONT33 \
    --ik1 ../../../Illumina/HepG2/Kallisto/Rep1/abundance.tsv \
    --ik2 ../../../Illumina/HepG2/Kallisto/Rep2/abundance.tsv \
    --dtype ONT \
    --color green \
    --outdir ../HepG2/
```
