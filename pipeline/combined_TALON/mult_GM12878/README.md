## Run our four GM12878 PacBio datasets through TALON in order to track gene detection by read depth.
Datasets:
- D6
- D7
- D8 (GM12878 Rep 1)
- D9 (GM12878 Rep 2)

Start from Tier 1 database.
```
cp ../full_gencode_v29_2019-03-12.db mult_GM12878.db
```

Then run TALON. 
```
./run_talon.sh
```

## Create unfiltered abundance file
```
python /pub/dwyman/TALON/post-TALON_tools/create_abundance_file_from_database.py \
           --db mult_GM12878.db \
           -a gencode_v29 \
           --build hg38 \
           --o mult_GM12878
```

## Plot gene detection as a function of read depth
```
Rscript ../../../analysis_scripts/plot_discovery_curve_knownOnly.R --f mult_GM12878_talon_abundance.tsv --color blue --rc read_counts.csv --ik1 ../../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv --ik2 ../../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv -o .
```
