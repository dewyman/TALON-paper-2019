# Post-TALON analysis on GM12878

## Make a whitelist file of filtered HepG2 transcripts
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
          --db full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist GM12878_whitelist.csv \
          --datasets GM12878_datasets.txt \
          --o GM12878_filtered
```

## Make a filtered abundance matrix
```
python /pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db full_gencode_v29_2019-03-12.db \
           -a gencode_v29 \
           --filter \
           -p pairings.csv \
           --o GM12878
```

## Make an abundance matrix without filtering
```
python /pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db full_gencode_v29_2019-03-12.db \
           -a gencode_v29 \
           --o GM12878
```

## Detection by TPM (without whitelist)
```
module load R/3.5.1
Rscript ../../../analysis_scripts/plot_detection_by_TPM_for_datasets.R \
           --db full_gencode_v29_2019-03-12.db \
           --datasets D8,D9 \
           --ik ../../../Illumina/GM12878/Kallisto/abundance.tsv \
           --color blue \
           -o GM12878_plots
```
