# Post-TALON analysis on K562

## Summary
```
python /pub/dwyman/TALON/post-TALON_tools/summarize_datasets.py \
          --db full_gencode_v29_2019-03-12.db \
          --v \
          --o K562
```

## Make a whitelist file of filtered K562 transcripts
```
python /dfs2/pub/dwyman/TALON/post-TALON_tools/filter_talon_transcripts.py \
          --db full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -p pairings.csv \
          --o K562_whitelist.csv
```

## Make a GTF file using the K562 whitelist
```
python /dfs2/pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist K562_whitelist.csv \
          --datasets K562_datasets.txt \
          --o K562_filtered
```

## Make a filtered abundance matrix  
```
python /dfs2/pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db full_gencode_v29_2019-03-12.db \
           -a gencode_v29 \
           --filter \
           -p pairings.csv \
           --o K562 
```
## Unfiltered abundance matrix
```
python /dfs2/pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db full_gencode_v29_2019-03-12.db \
           -a gencode_v29 \
           --o K562
```


## Detection by TPM (without whitelist)
```
module load R/3.5.1
Rscript ../../../analysis_scripts/plot_detection_by_TPM_for_datasets.R \
           --db full_gencode_v29_2019-03-12.db \
           --datasets D10,D11 \
           --ik ../../../Illumina/K562/Kallisto/abundance.tsv \
           --color red \
           -o K562_plots
```

## Plot gene correlation of D10 and D11
```
module load R/3.5.1
Rscript ../../../analysis_scripts/plot_pacbio_gene_expression_corr.R \
          --f K562_talon_abundance.tsv \
          --color blue \
          --d1 D10 --d2 D11 \
          --celltype K562 \
          -o K562_plots

Rscript ../../../analysis_scripts/plot_pacbio_gene_expression_corr.R \
          --f K562_talon_abundance.tsv \
          --color blue \
          --d1 D10 --d2 D11 \
          --celltype K562 \
          --intergenic \
          -o K562_plots

Rscript ../../../analysis_scripts/plot_pacbio_gene_expression_corr.R \
          --f K562_talon_abundance.tsv \
          --color blue \
          --d1 D10 --d2 D11 \
          --celltype K562 \
          --antisense \
          --intergenic \
          -o K562_plots
```

## Plot transcript correlation
```
Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f K562_talon_abundance_filtered.tsv \
          --d1 D10 --d2 D11 \
          --celltype K562 \
          -o K562_plots

Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f K562_talon_abundance_filtered.tsv \
          --d1 D10 --d2 D11 \
          --celltype K562 \
          --ISM \
          -o K562_plots

Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f K562_talon_abundance_filtered.tsv \
          --d1 D10 --d2 D11 \
          --celltype K562 \
          --NIC \
          -o K562_plots

Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f K562_talon_abundance_filtered.tsv \
          --d1 D10 --d2 D11 \
          --celltype K562 \
          --NNC \
          -o K562_plots

Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f K562_talon_abundance_filtered.tsv \
          --d1 D10 --d2 D11 \
          --celltype K562 \
          --ISM \
          --NIC \
          --NNC \
          --antisense \
          --intergenic \
          --genomic \
          -o K562_plots
```

## MA plot for known genes
```
Rscript /pub/dwyman/TALON-paper-2019/analysis_scripts/plot_TPM_chisquare_pvalues.R \
    --f K562_talon_abundance.tsv \
    --datasets D10,D11 \
    --ik ../../../Illumina/K562/Kallisto/abundance.tsv \
    --color green \
    -o K562_plots
```

## MA plot for known transcripts
```
Rscript /pub/dwyman/TALON-paper-2019/analysis_scripts/MA_plot_for_transcripts.R \
    --f K562_talon_abundance.tsv \
    --datasets D10,D11 \
    --ik ../../../Illumina/K562/Kallisto/abundance.tsv \
    --color green \
    -o K562_plots
```


## Novelty categories plots
```
Rscript ../../../analysis_scripts/plot_novelty_categories.R \
        --db full_gencode_v29_2019-03-12.db \
        --w K562_whitelist.csv \
        --datasets D10,D11 \
        -o K562_plots
```

## RNA-PET analysis
```
mkdir -p RNA-PET
module load bedtools
cd ../../../RNA-PET
python run_RNA-PET_analysis.py \
    --gtf ../pipeline/combined_TALON/K562/K562_filtered_talon.gtf \
    --rnapet data/K562_hg38.bed \
    --maxdist 100 \
    --o ../pipeline/combined_TALON/K562/RNA-PET/K562

Rscript plot_RNA-PET_support.R \
    --f ../pipeline/combined_TALON/K562/RNA-PET/K562_RNA-PET_results.csv \
    --novelty ../pipeline/combined_TALON/K562/RNA-PET/transcript_beds/K562_novelty.csv \
    -o ../pipeline/combined_TALON/K562/RNA-PET/K562
```

## Compare long read K562 splice jns to K562 short reads
```
module load dwyman/anaconda/3
source activate runTC
python /data/users/dwyman/TranscriptClean-1.0.7/accessory_scripts/get_SJs_from_gtf.py \
    --f K562_filtered_talon.gtf \
    --g ../../../refs/hg38/hg38.fa \
    --o K562_filtered_talon_spliceJns.txt
source deactivate

python ../../../analysis_scripts/compare_sjs/compare_sjs.py \
    --short ../../../Illumina/K562/STAR/SJ.out.tab \
    --long K562_filtered_talon_spliceJns.txt \
    --o K562_filtered_talon
```
