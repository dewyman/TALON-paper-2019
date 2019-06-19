# set -e 
bash generate_abundance_tables.sh
bash generate_gtfs.sh
bash generate_plots.hs
bash rename_abundance_file_datasets.sh