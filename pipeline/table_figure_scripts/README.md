Update local bash variables to appropriate data/TALON/TALON-paper-2019 directory in each script before running!

## Generate abundance tables (used for plotting and in the supplement)
```
bash generate_abundance_tables.sh
```

## Generate GTFs (used for track generation and in the supplement)
```
bash generate_gtfs.sh
```

## Generate plots seen in the main paper and in the supplement
```
bash generate_plots.sh
```

## Rename abundance file dataset name columns to be more descriptive
```
bash rename_abundance_file_datasets.sh
```

## TALONclass.py
```
Used by scripts S35table.py and S36table.py
```

## S35table.py
```
Generates a table with gene read counts based on novelty category
```

## S36table.py
```
Generate table with genes with higher novelty counts for cortex and hippocampus, also generates a plot to show how many of them are unique or shared between tissues
```
