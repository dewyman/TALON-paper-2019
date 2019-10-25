# EBV

We wanted to see if TALON and PacBio long-read sequencing could be used to detect, characterize, and quantify transcripts from the EBV chromosome used to immortalize the GM12878 cell line. 


## TALON
### Initialize TALON db with custom EBV gtf (using chr1 because sam files are fake mapped to chr1)
```bash
sed 's/chrEBV/chr1/' ebv.gtf > ebv_chr1.gtf
talon_initialize_database \
    --f ebv_chr1.gtf \
    --g HHV4 \
    --a HHV4 \
    --o ebv
```

### Run TALON on EBV reads extracted from mapped GM12878 sam files
```bash
talon \
    --f ebv_talon_config.csv \
    --db ebv.db \
    --build HHV4 \
    --o ebv
```

## Post-TALON files

### Get a whitelist of transcripts via TALON filtering
```bash
printf "PacBio_GM12878_1,PacBio_GM12878_2" > pairings
talon_filter_transcripts \
    --db ebv.db \
    -a HHV4 \
    -p pairings \
    --o ebv_whitelist
```

### Get filtered GTF file
```bash
talon_create_GTF \
      --db ebv.db \
      --b HHV4 \
      --a HHV4 \
      --o ebv \
      --whitelist ebv_whitelist \
      --observed
```

### Get unfiltered abundance file
```bash
talon_abundance \
    --db ebv.db \
    --a HHV4 \
    --b HHV4 \
    --o ebv
```

### Get filtered abundance file
```bash
talon_abundance \
    --db ebv.db \
    --a HHV4 \
    --b HHV4 \
    --whitelist ebv_whitelist \
    --o ebv
```

### GM12878 files to compare with 
```bash
sup_tables=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_10-19/human_TALON/analysis/supplementary_tables/
ln -s ${sup_tables}S2_GM12878_talon_observedOnly.gtf S2_GM12878_talon_observedOnly.gtf
ln -s ${sup_tables}S3_GM12878_talon_abundance.tsv S3_GM12878_talon_abundance.tsv
ln -s ${sup_tables}S4_GM12878_talon_abundance_filtered.tsv S4_GM12878_talon_abundance_filtered.tsv
```

## Plotting
### Generate EBV abundance violin plots
```bash
python ebv_compute_tpms.py --c ebv_expression_config.csv
Rscript plot_ebv_v_human_abundances.R \
          --gene_csv ebv_human_gene_abundance.csv \
          --transcript_csv ebv_human_transcript_abundance.csv \
          --datasets combined
```

<img align="center" width="300" src="combined_genes_ebv_human.png "><img align="center" width="300" src="combined_transcripts_ebv_human.png ">

## Genome browser trackline
### Generate tracklines using above GTF
```bash
python ../analysis_scripts/gen_novelty_tracks_gtf.py \
          --c ebv_gtf_track_config.csv
url=`cut -d, -f5 ebv_gtf_track_config.csv`
n=`cut -d, -f2 ebv_gtf_track_config.csv`
cp ebv_chr1.gtf ebv_talon_observedOnly_tracks/
printf 'track name="EBV Reference" visibility=pack color=0,0,128\n%s/ebv_chr1.gtf' "$url" >> ebv_talon_observedOnly_tracks/ebv_talon_observedOnly_${n}_tracks
```

From here, you can open the genome browser up and display your tracklines, and use the genome browser's PDF functionality or take a screenshot to get the genome browser shot. 

<img align="center" width="700" src="ebv_browser.png ">

