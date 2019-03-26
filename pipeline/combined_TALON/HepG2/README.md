# Post-TALON analysis on HepG2

## Summary
```
python /pub/dwyman/TALON/post-TALON_tools/summarize_datasets.py \
          --db ../full_gencode_v29_2019-03-12.db \
          --v \
          --o HepG2
```

## Make a whitelist file of filtered HepG2 transcripts
```
python /pub/dwyman/TALON/post-TALON_tools/filter_talon_transcripts.py \
          --db ../full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -p pairings.csv \
          --o HepG2_whitelist.csv
```

## Make a GTF file using the HepG2 whitelist
```
python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist HepG2_whitelist.csv \
          --datasets HepG2_datasets.txt \
          --o HepG2_filtered
```

## Make a filtered abundance matrix  
```
python /pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db ../full_gencode_v29_2019-03-12.db \
           -a gencode_v29 \
           --filter \
           -p pairings.csv \
           --o HepG2 
```

## Detection by TPM (without whitelist)
```
module load R/3.5.1
Rscript ../../../analysis_scripts/plot_detection_by_TPM_for_datasets.R \
           --db ../full_gencode_v29_2019-03-12.db \
           --datasets D4,D5 \
           --ik ../../../Illumina/HepG2/Kallisto/abundance.tsv \
           --color green \
           -o HepG2_plots
```

## Plot correlation of D4 and D5
```
module load R/3.5.1
Rscript ../../../analysis_scripts/plot_pacbio_expression_corr.R \
          --f HepG2_talon_abundance_filtered.tsv \
          --w HepG2_whitelist.csv \
          --color green \
          --d1 D4 --d2 D5 \
          -o HepG2_plots

Rscript ../../../analysis_scripts/plot_pacbio_gene_expression_corr.R \
          --f HepG2_talon_abundance.tsv \
          --color green \
          --d1 D4 --d2 D5 \
          -o HepG2_plots

Rscript ../../../analysis_scripts/plot_pacbio_transcript_expression_corr.R \
          --f HepG2_talon_abundance.tsv \
          --d1 D4 --d2 D5 \
          -o HepG2_plots
```

## MA plot for known genes
```
Rscript /pub/dwyman/TALON-paper-2019/analysis_scripts/plot_TPM_chisquare_pvalues.R \
    --f HepG2_talon_abundance.tsv \
    --datasets D4,D5 \
    --ik ../../../Illumina/HepG2/Kallisto/abundance.tsv \
    --color green \
    -o HepG2_plots
```

## MA plot for known transcripts
```
Rscript /pub/dwyman/TALON-paper-2019/analysis_scripts/MA_plot_for_transcripts.R \
    --f HepG2_talon_abundance.tsv \
    --datasets D4,D5 \
    --ik ../../../Illumina/HepG2/Kallisto/abundance.tsv \
    --color green \
    -o HepG2_plots
```


## Novelty categories plots
```
Rscript /pub/dwyman/TALON-paper-2019/analysis_scripts/plot_novelty_categories.R \
        --db ../full_gencode_v29_2019-03-12.db \
        --w HepG2_whitelist.csv \
        -o HepG2_plots
```
