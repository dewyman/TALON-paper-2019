# Post-TALON analysis on HepG2

## Make a whitelist file of filtered HepG2 transcripts
```
python ~/TALON-4.0/post-TALON_tools/filter_talon_transcripts.py \
          --db ../full_gencode_v29_2019-03-12.db \
          -a gencode_v29 \
          -p pairings.csv
          --o HepG2_whitelist.csv
```
