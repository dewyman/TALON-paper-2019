
```
cp /pub/dwyman/TALON-paper-2019/refs/TALON/unmodified_full_gencode_v29_2019-03-12.db full_gencode_v29_2019-03-12.db
```

After running TALON on all 3 cell lines:

## Make a whitelist file of transcripts for all 3 cell lines
```
python /dfs2/pub/dwyman/TALON/post-TALON_tools/filter_talon_transcripts.py \
          --db full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -p pairings.csv \
          --o tier1_whitelist.csv
```

## Make a GTF file using the GM12878 whitelist
```
python /dfs2/pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist tier1_whitelist.csv \
          --datasets tier1_datasets.txt \
          --o tier1_filtered
```

