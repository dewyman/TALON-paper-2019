# EBV

We wanted to see if TALON and PacBio long-read sequencing could be used to detect, characterize, and quantify transcripts from the EBV chromosome used to immortalize the GM12878 cell line. 


## TALON
### Initialize TALON db with custom EBV gtf (using chr1 because sam files are fake mapped to chr1)
```
sed 's/chrEBV/chr1/' ebv.gtf > ebv_chr1.gtf
python ~/mortazavi_lab/bin/TALON/initialize_talon_database.py \
          --f ebv_chr1.gtf \
          --g HHV4 \
          --a HHV4 \
          --o ebv
```

### Run TALON on EBV reads extracted from mapped GM12878 sam files
```
python ~/mortazavi_lab/bin/TALON/talon.py \
          --f ebv_talon_config.csv \
          --db ebv.db \
          --build HHV4 \
          --o ebv
```

## Post-TALON files
### Get unfiltered abundance file
```
python ~/mortazavi_lab/bin/TALON/post-TALON_tools/create_abundance_file_from_database.py \
          --db ebv.db \
          --a HHV4 \
          --b HHV4 \
          --o ebv
```

### Get filtered abundance file
```
printf "D8,D9" > pairings
python ~/mortazavi_lab/bin/TALON/post-TALON_tools/create_abundance_file_from_database.py \
          --db ebv.db \
          --a HHV4 \
          --b HHV4 \
          --filter \
          --p pairings \
          --o ebv
```

### Get unfiltered whole GM12878 abundance file
```
printf "D8\nD9" > datasets
python ~/mortazavi_lab/bin/TALON/post-TALON_tools/create_abundance_file_from_database.py \
          --db ~/mortazavi_lab/data/ont_tier1/full_gencode_v29_2019-05-24.db \
          --a gencode_v29 \
          --b hg38 \
          --o gm12878
python ../../analysis_scripts/get_dataset_specific_abundance.py \
          --infile gm12878_talon_abundance.tsv \
          --d datasets
```

### Get filtered whole GM12878 abundance file
```
python ~/mortazavi_lab/bin/TALON/post-TALON_tools/create_abundance_file_from_database.py \
          --db ~/mortazavi_lab/data/ont_tier1/full_gencode_v29_2019-05-24.db \
          --a gencode_v29 \
          --b hg38 \
          --o gm12878 \
          --filter \
          --pairings pairings
python ../../analysis_scripts/get_dataset_specific_abundance.py \
          --infile gm12878_talon_abundance_filtered.tsv \
          --d datasets
```

### Get GTF file
```
python ~/mortazavi_lab/bin/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ebv.db \
          --b HHV4 \
          --a HHV4 \
          --o ebv \
          --observed
```

### Get GM12878 GTF file
```
python ~/mortazavi_lab/bin/TALON/post-TALON_tools/create_GTF_from_database.py \
            --db ~/mortazavi_lab/data/ont_tier1/full_gencode_v29_2019-05-24.db \
            --a gencode_v29 \
            --b hg38 \
            --o gm12878 \
            --observed
```

## Genome browser trackline
### Generate tracklines using above GTF
```
python ../../analysis_scripts/gen_novelty_tracks_gtf.py \
          --c ebv_gtf_track_config.csv
url=`cut -d, -f5 ebv_gtf_track_config.csv`
n=`cut -d, -f2 ebv_gtf_track_config.csv`
cp ebv_chr1.gtf ebv_talon_observedOnly_tracks/
printf 'track name="EBV Reference" visibility=pack color=0,0,128\n%s/ebv_chr1.gtf' "$url" >> ebv_talon_observedOnly_tracks/ebv_talon_observedOnly_${n}_tracks
```

## Plotting
### Generate EBV abundance violin plots
```
python ebv_compute_tpms.py --c ebv_expression_config.csv
Rscript plot_ebv_v_human_abundances.R \
          --gene_csv ebv_human_gene_abundance.csv \
          --transcript_csv ebv_human_transcript_abundance.csv\
          --datasets combined
```
