import pandas as pd 
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description=\
	'Renames PacBio and ONT datasets with more\
	 intelligent names')
parser.add_argument('--f', help='file to swap dataset col names in')
args = parser.parse_args()
f = args.f

# read in mapping file
map_df = pd.read_csv('dataset_id_name_map.tsv', sep='\t')
map_df.set_index('dataset_id', inplace=True)
map_dict = map_df.to_dict()

df = pd.read_csv(f, sep='\t')
df.rename(columns=map_dict['dataset_name'], inplace=True)
df.to_csv(f, sep='\t', index=False)


