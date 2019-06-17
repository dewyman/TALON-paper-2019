import pandas as pd 
from collections import defaultdict
import glob

# read in mapping file
map_df = pd.read_csv('dataset_id_name_map.tsv', sep='\t')
map_df.set_index('dataset_id', inplace=True)
map_dict = map_df.to_dict()

for f in glob.glob('*talon_abundance*tsv'):
	df = pd.read_csv(f, sep='\t')
	df.rename(columns=map_dict['dataset_name'], inplace=True)
	df.to_csv(f, sep='\t', index=False)


