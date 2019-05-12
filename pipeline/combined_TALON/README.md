
```
cp /pub/dwyman/TALON-paper-2019/refs/TALON/unmodified_full_gencode_v29_2019-03-12.db full_gencode_v29_2019-03-12.db
```

After running TALON on all 3 cell lines:

## Summary
```
python /pub/dwyman/TALON/post-TALON_tools/summarize_datasets.py \
          --db full_gencode_v29_2019-03-12.db \
          --groups summary_groups.csv \
          --v \
          --o tier1
```

## Make a whitelist file of transcripts for all 3 cell lines, allowing different cell liens to corroborate each other
```
python /pub/dwyman/TALON/post-TALON_tools/filter_talon_transcripts.py \
          --db full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -p cross-cell-pairings.csv \
          --o tier1_whitelist.csv
```

## Make a GTF file using the cross-cell-line whitelist
```
python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist tier1_whitelist.csv \
          --datasets tier1_datasets.txt \
          --o tier1_filtered
```

## Make color-coded GTF tracks
```
python ../../analysis_scripts/gen_novelty_tracks_gtf.py --c track_config.csv
```

## Make abundance matrix
```
python /pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db full_gencode_v29_2019-03-12.db \
           -a gencode_v29 \
           --build hg38 \
           --o tier1
```

## Make filtered abundance matrix
```
python /pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db full_gencode_v29_2019-03-12.db \
           -a gencode_v29 \
           --filter \
           --build hg38 \
           -p cross-cell-pairings.csv \
           --o tier1
```
## Make a fasta file from the filtered GTF
```
module load perl
module load freese/TransDecoder
export PERL_HASH_SEED=0
UTILPATH=/data/apps/user_contributed_software/freese/TransDecoder/5.5.0/util
$UTILPATH/gtf_genome_to_cdna_fasta.pl \
    tier1_filtered_talon.gtf \
    ../../refs/hg38/hg38.fa > tier1_filtered_talon.fasta
```


## Make abundance plots for TCF3 transcripts (for genome browser plot)
```
Rscript ../../analysis_scripts/plot_expression_for_genome_browser.R \
        --f tier1_talon_abundance_filtered.tsv \ 
        --groups groups.csv \
        --transcripts TCF3_transcript_names.txt \
        -o tier1_plots
```
