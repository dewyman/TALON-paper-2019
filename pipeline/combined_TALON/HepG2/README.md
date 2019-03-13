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
python /pub/dwyman/post-TALON_tools/filter_talon_transcripts.py \
          --db ../full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -p pairings.csv \
          --o HepG2_whitelist.csv
```

## Make a GTF file using the HepG2 whitelist
```
python /pub/dwyman/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist HepG2_whitelist.csv \
          --datasets HepG2_datasets.txt \
          --o HepG2_filtered
```

## Make a filtered abundance matrix  
```
python /pub/dwyman/post-TALON_tools/create_abundance_file_from_database.py \
           --db ../full_gencode_v29_2019-03-12.db \
           -a gencode_v29 \
           --filter \
           --o HepG2 
```
