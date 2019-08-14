## Make a whitelist file of filtered GM12878 transcripts
```
python /pub/dwyman/TALON/post-TALON_tools/filter_talon_transcripts.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -p pairings.csv \
          --o GM12878_whitelist.csv
```

## Make a GTF file using the GM12878 whitelist
```
python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist GM12878_whitelist.csv \
          --datasets GM12878_datasets.txt \
          --o GM12878_filtered
```

## Make an abundance matrix without filtering
```
python /pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db ../full_gencode_v29_2019-05-24.db \
           -a gencode_v29 \
           --build hg38 \
           --o ../all
```

## Make a filtered abundance matrix 
```
python ~/mortazavi_lab/bin/TALON/post-TALON_tools/create_abundance_file_from_database.py \
    --db ~/mortazavi_lab/data/ont_tier1/full_gencode_v29_2019-05-24.db \
    -a gencode_v29 \
    --build hg38 \
    --o all \
    --filter \
    -p pairings.csv \
    --o all
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
    --gtf ../pipeline/nanopore/GM12878/GM12878_filtered_talon.gtf \
    --rnapet data/GM12878_hg38.bed \
    --maxdist 100 \
    --o ../pipeline/nanopore/GM12878/RNA-PET/GM12878

cd ../pipeline/nanopore/GM12878
Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R  \
     --f RNA-PET/GM12878_RNA-PET_results.csv   \
     --t RNA-PET \
     --novelty RNA-PET/transcript_beds/GM12878_novelty.csv \
     --abundance ../all_talon_abundance.tsv \
     --d1 ONT33 --d2 ONT34 \
     --as ../all_antisense_mapping.csv \
     --splitISM \
     -o RNA-PET/GM12878
```

## CAGE analysis, ENCODE
```
source activate mypython3.7.2
mkdir -p CAGE/ENCODE/
python ../../../CAGE/run_CAGE_analysis.py \
        --gtf GM12878_filtered_talon.gtf \
        --cage ../../../CAGE/data/ENCODE/GM12878_CAGE.bed \
        --maxdist 100 \
        --o CAGE/ENCODE/GM12878

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R  \
     --f CAGE/ENCODE/GM12878_CAGE_results.csv   \
     --t CAGE \
     --novelty CAGE/ENCODE/transcript_beds/GM12878_novelty.csv \
     --abundance ../all_talon_abundance.tsv \
     --d1 ONT33 --d2 ONT34 \
     --as ../all_antisense_mapping.csv \
     --splitISM \
     -o CAGE/ENCODE/GM12878
```

## CAGE analysis, FANTOM
```
source activate mypython3.7.2
mkdir -p CAGE/FANTOM5/
python ../../../CAGE/run_CAGE_analysis.py \
        --gtf GM12878_filtered_talon.gtf \
        --cage ../../../CAGE/data/FANTOM5/hg38_CAGE.bed \
        --maxdist 100 \
        --o CAGE/FANTOM5/GM12878

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R  \
     --f CAGE/FANTOM5/GM12878_CAGE_results.csv   \
     --t CAGE \
     --novelty CAGE/FANTOM5/transcript_beds/GM12878_novelty.csv \
     --abundance ../all_talon_abundance.tsv \
     --d1 ONT33 --d2 ONT34 \
     --as ../all_antisense_mapping.csv \
     --splitISM \
     -o CAGE/FANTOM5/GM12878
```

## Computational PAS analysis
```
source activate mypython3.7.2
mkdir PAS-comp
python ../../../PAS-computational/run_computational_PAS_analysis.py \
        --gtf GM12878_filtered_talon.gtf \
        --genome ../../../refs/hg38/hg38.fa \
        --maxdist 35 \
        --o PAS-comp/GM12878

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R \
    --f PAS-comp/GM12878_polyA_motif.csv \
    --t PAS-comp \
    --novelty PAS-comp/transcript_beds/GM12878_novelty.csv \
    --abundance ../all_talon_abundance.tsv \
    --d1 ONT33 --d2 ONT34 --as ../all_antisense_mapping.csv \
    --splitISM \
    -o PAS-comp/GM12878
```

## ONT vs. Illumina quantification analysis
```
# gene detection
Rscript ../../../compare_to_FLAIR/plot_detection_by_TPM_for_datasets.R \
    --f all_talon_abundance.tsv \
    --datasets ONT24,ONT25 \
    --ik1 ../../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
    --ik2 ../../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
    --color blue \
    --outdir ../GM12878/

# differences in gene quantification
Rscript ../../../analysis_scripts/pacbio_v_illumina_edgeR.R \
    --f all_talon_abundance.tsv \
    --datasets ONT24,ONT25 \
    --ik1 ../../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
    --ik2 ../../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
    --dtype ONT \
    --color blue \
    --outdir ../GM12878/

# differences in transcript quantification 
Rscript ../../../analysis_scripts/pacbio_v_illumina_edgeR_transcripts.R \
    --f all_talon_abundance.tsv \
    --datasets ONT24,ONT25 \
    --ik1 ../../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
    --ik2 ../../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
    --dtype ONT \
    --color blue \
    --outdir ../GM12878/
```
