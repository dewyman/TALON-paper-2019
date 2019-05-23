# Post-TALON analysis on GM12878

## Make a whitelist file of filtered GM12878 transcripts
```
python /pub/dwyman/TALON/post-TALON_tools/filter_talon_transcripts.py \
          --db full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -p pairings.csv \
          --o GM12878_whitelist.csv
```

## Make a GTF file using the GM12878 whitelist
```
python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist GM12878_whitelist.csv \
          --datasets GM12878_datasets.txt \
          --o GM12878_filtered
```

## Make a filtered abundance matrix
```
python /pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db ../full_gencode_v29_2019-03-12.db \
           -a gencode_v29 \
           --build hg38 \
           --filter \
           -p pairings.csv \
           --o GM12878
```

## Make an abundance matrix without filtering
```
python /pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db ../full_gencode_v29_2019-03-12.db \
           -a gencode_v29 \
           --build hg38 \
           --o GM12878
```

## Detection by TPM (without whitelist)
```
module load R/3.5.1
Rscript ../../../analysis_scripts/plot_detection_by_TPM_for_datasets.R \
           --db ../full_gencode_v29_2019-03-12.db \
           --datasets D8,D9 \
           --ik1 ../../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
           --ik2 ../../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
           --color blue \
           -o GM12878_plots
```
## Gene length by detection
```
Rscript ../../../analysis_scripts/plot_gene_length_by_detection_for_datasets.R \
           --db ../full_gencode_v29_2019-03-12.db \
           --datasets D8,D9 \
           --ik1 ../../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
           --ik2 ../../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
           --color blue \
           -o GM12878_plots
```

## Plot gene correlation of D8 and D9
```
module load R/3.5.1
Rscript ../../../analysis_scripts/plot_pacbio_gene_expression_corr.R \
          --f GM12878_talon_abundance.tsv \
          --color blue \
          --d1 D8 --d2 D9 \
          --celltype GM12878 \
          -o GM12878_plots

Rscript ../../../analysis_scripts/plot_pacbio_gene_expression_corr.R \
          --f GM12878_talon_abundance.tsv \
          --color blue \
          --d1 D8 --d2 D9 \
          --celltype GM12878 \
          --intergenic \
          -o GM12878_plots

Rscript ../../../analysis_scripts/plot_pacbio_gene_expression_corr.R \
          --f GM12878_talon_abundance.tsv \
          --color blue \
          --d1 D8 --d2 D9 \
          --celltype GM12878 \
          --intergenic \
          --antisense \
          -o GM12878_plots

Rscript ../../../analysis_scripts/plot_pacbio_gene_expression_corr.R \
          --f GM12878_talon_abundance.tsv \
          --color blue \
          --d1 D8 --d2 D9 \
          --celltype GM12878 \
          --antisense \
          -o GM12878_plots

```

## Plot transcript correlation
```
Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f GM12878_talon_abundance_filtered.tsv \
          --d1 D8 --d2 D9 \
          --celltype GM12878 \
          -o GM12878_plots

Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f GM12878_talon_abundance_filtered.tsv \
          --d1 D8 --d2 D9 \
          --celltype GM12878 \
          --ISM \
          -o GM12878_plots

Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f GM12878_talon_abundance_filtered.tsv \
          --d1 D8 --d2 D9 \
          --celltype GM12878 \
          --NIC \
          -o GM12878_plots

Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f GM12878_talon_abundance_filtered.tsv \
          --d1 D8 --d2 D9 \
          --celltype GM12878 \
          --NNC \
          -o GM12878_plots

Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f GM12878_talon_abundance_filtered.tsv \
          --d1 D8 --d2 D9 \
          --celltype GM12878 \
          --ISM \
          --NIC \
          --NNC \
          -o GM12878_plots

Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f GM12878_talon_abundance_filtered.tsv \
          --d1 D8 --d2 D9 \
          --celltype GM12878 \
          --NIC \
          --NNC \
          -o GM12878_plots

Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f GM12878_talon_abundance_filtered.tsv \
          --d1 D8 --d2 D9 \
          --celltype GM12878 \
          --ISM \
          --NIC \
          --NNC \
          --antisense \
          --intergenic \
          -o GM12878_plots
```
## Demo scatterplot just for Figure 1A
```
Rscript ../../../analysis_scripts/plot_simplified_corr_fig1.R \
          --f GM12878_talon_abundance.tsv \
          --d1 D8 --d2 D9 \
          --color blue \
          -o GM12878_plots
```

## edgeR on Illumina vs PacBio, two bioreps each
```
Rscript ../../../analysis_scripts/pacbio_v_illumina_edgeR.R \
    --f GM12878_talon_abundance.tsv \
    --datasets D8,D9 \
    --ik1 ../../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
    --ik2 ../../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
    --color blue \
    -o GM12878_plots 

Rscript ../../../analysis_scripts/pacbio_v_illumina_edgeR_transcripts.R \
    --f GM12878_talon_abundance.tsv \
    --datasets D8,D9 \
    --ik1 ../../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
    --ik2 ../../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
    --color green \
    -o GM12878_plots
```

## Novelty categories plots
```
Rscript ../../../analysis_scripts/plot_novelty_categories.R \
        --db full_gencode_v29_2019-03-12.db \
        --w GM12878_whitelist.csv \
        --datasets D8 \
        -o GM12878_plots
```

## Map antisense TALON genes to their sense counterparts
```
python /dfs2/pub/dwyman/TALON/post-TALON_tools/map_antisense_genes_to_sense.py \
    --db full_gencode_v29_2019-03-12.db \
    --a gencode_v29 
    --o GM12878
```

## RNA-PET analysis
```
mkdir -p RNA-PET
module load bedtools
cd ../../../RNA-PET
python run_RNA-PET_analysis.py \
    --gtf ../pipeline/combined_TALON/GM12878/GM12878_filtered_talon.gtf \
    --rnapet data/GM12878_hg38.bed \
    --maxdist 100 \
    --o ../pipeline/combined_TALON/GM12878/RNA-PET/GM12878
cd ../pipeline/combined_TALON/GM12878
Rscript ../../../RNA-PET/plot_RNA-PET_support.R  \
     --f RNA-PET/GM12878_RNA-PET_results.csv   \
     --novelty RNA-PET/transcript_beds/GM12878_novelty.csv \
     --abundance GM12878_talon_abundance.tsv \
     --d1 D8 --d2 D9 \
     --as GM12878_antisense_mapping.csv \
     -o RNA-PET/GM12878
```

## CAGE analysis, ENCODE
```
source activate mypython3.7.2
python ../../../CAGE/run_CAGE_analysis.py \
        --gtf GM12878_filtered_talon.gtf \
        --cage ../../../CAGE/data/ENCODE/GM12878_CAGE.bed \
        --maxdist 100 \
        --o CAGE/ENCODE/GM12878

Rscript ../../../CAGE/plot_CAGE_support.R \
    --f CAGE/ENCODE/GM12878_CAGE_results.csv \
    --novelty CAGE/ENCODE/transcript_beds/GM12878_novelty.csv \
    --abundance GM12878_talon_abundance.tsv \
    --d1 D8 --d2 D9 --as GM12878_antisense_mapping.csv \
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

Rscript ../../../CAGE/plot_CAGE_support.R \
    --f CAGE/FANTOM5/GM12878_CAGE_results.csv \
    --novelty CAGE/FANTOM5/transcript_beds/GM12878_novelty.csv \
    --abundance GM12878_talon_abundance.tsv \
    --d1 D8 --d2 D9 --as GM12878_antisense_mapping.csv \
    -o CAGE/FANTOM5/GM12878
```
## PAS analysis, GENCODE manual PolyA annotation
```
source activate mypython3.7.2
python ../../../PAS-seq/run_GENCODE_PAS-seq_analysis.py \
        --gtf GM12878_filtered_talon.gtf \
        --pas ../../../PAS-seq/gencode.v29.metadata.PolyA_feature.bed \
        --maxdist 35 \
        --o PAS-annot/GM12878
```

## Compare long read GM12878 splice jns to GM12878 short reads
```
module load dwyman/anaconda/3
source activate runTC
python /data/users/dwyman/TranscriptClean-1.0.7/accessory_scripts/get_SJs_from_gtf.py \
    --f GM12878_filtered_talon.gtf \
    --g ../../../refs/hg38/hg38.fa \
    --o GM12878_filtered_talon_spliceJns.txt
source deactivate

python ../../../analysis_scripts/compare_sjs/compare_sjs.py \
    --short ../../../Illumina/GM12878/STAR/SJ.out.tab \
    --long GM12878_filtered_talon_spliceJns.txt \
    --o GM12878_filtered_talon
```
