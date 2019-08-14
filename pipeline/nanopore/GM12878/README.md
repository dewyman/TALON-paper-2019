## Make a whitelist file of filtered GM12878 transcripts
```
python /pub/dwyman/TALON/post-TALON_tools/filter_talon_transcripts.py \
          --db ../hg38_v29_incomplete.db \
          -a v29 \
          -p pairings.csv \
          --o GM12878_whitelist.csv
```

## Make a GTF file using the GM12878 whitelist
```
python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../hg38_v29_incomplete.db \
          -a v29 \
          -b hg38 \
          --whitelist GM12878_whitelist.csv \
          --datasets GM12878_datasets.txt \
          --o GM12878_filtered
```

## Make an abundance matrix without filtering
```
python /pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db ../hg38_v29_incomplete.db \
           -a v29 \
           --build hg38 \
           --o GM12878
```

## Map antisense TALON genes to their sense counterparts
```
python /pub/dwyman/TALON/post-TALON_tools/map_antisense_genes_to_sense.py \
    --db ../hg38_v29_incomplete.db \
    --a v29 \
    --o GM12878
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
    --o ../pipeline/nanopore/RNA-PET/GM12878/GM12878

cd ../pipeline/nanopore/GM12878
Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R  \
     --f RNA-PET/GM12878_RNA-PET_results.csv   \
     --t RNA-PET \
     --novelty RNA-PET/transcript_beds/GM12878_novelty.csv \
     --abundance GM12878_talon_abundance.tsv \
     --d1 ONT24 --d2 ONT25 \
     --as GM12878_antisense_mapping.csv \
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
     --abundance GM12878_talon_abundance.tsv \
     --d1 ONT24 --d2 ONT25 \
     --as GM12878_antisense_mapping.csv \
     --splitISM \
     -o CAGE/ENCODE/GM12878
```

## CAGE analysis, FANTOM
```
source activate mypython3.7.2
python ../../../CAGE/run_CAGE_analysis.py \
        --gtf GM12878_filtered_talon.gtf \
        --cage ../../../CAGE/data/FANTOM5/hg38_CAGE.bed \
        --maxdist 100 \
        --o CAGE/FANTOM5/GM12878

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R  \
     --f CAGE/FANTOM5/GM12878_CAGE_results.csv   \
     --t CAGE \
     --novelty CAGE/FANTOM5/transcript_beds/GM12878_novelty.csv \
     --abundance GM12878_talon_abundance.tsv \
     --d1 ONT24 --d2 ONT25 \
     --as GM12878_antisense_mapping.csv \
     --splitISM \
     -o CAGE/FANTOM5/GM12878
```
## PAS analysis, GENCODE manual PolyA annotation
```
source activate mypython3.7.2
mkdir -p PAS-annot
python ../../../PAS-seq/run_GENCODE_PAS-seq_analysis.py \
        --gtf GM12878_filtered_talon.gtf \
        --pas ../../PAS-seq/gencode.v29.metadata.PolyA_feature.bed \
        --maxdist 50 \
        --o PAS-annot/GM12878

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R \
    --f PAS-annot/GM12878_PAS_results.csv \
    --t PAS-annot \
    --novelty PAS-annot/transcript_beds/GM12878_novelty.csv \
    --abundance GM12878_talon_abundance.tsv \
    --d1 ONT24 --d2 ONT25 --as GM12878_antisense_mapping.csv \
    --splitISM \
    -o PAS-annot/GM12878
```

## Computational PAS analysis
```
source activate mypython3.7.2
mkdir PAS-comp
python ../../../PAS-computational/run_computational_PAS_analysis.py \
        --gtf GM12878_filtered_talon.gtf \
        --genome ../../refs/hg38/hg38.fa \
        --maxdist 35 \
        --o PAS-comp/GM12878

Rscript ../../../analysis_scripts/plot_support_by_novelty_type.R \
    --f PAS-comp/GM12878_polyA_motif.csv \
    --t PAS-comp \
    --novelty PAS-comp/transcript_beds/GM12878_novelty.csv \
    --abundance GM12878_talon_abundance.tsv \
    --d1 ONT24 --d2 ONT25 --as GM12878_antisense_mapping.csv \
    --splitISM \
    -o PAS-comp/GM12878
```

