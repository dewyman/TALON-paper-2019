import pandas as pd 
import argparse

parser = argparse.ArgumentParser(description=\
		  'From an exisiting abundance file, cut out\
		   input dataset fields')
parser.add_argument('--f', help='abundance file')
parser.add_argument('--d', help='comma-separated string of desired dataset names')
args = parser.parse_args()

infile = args.f
datasets = args.d.split(',')

fields = ['gene_ID','transcript_ID','annot_gene_id',\
		  'annot_transcript_id', 'annot_gene_name', \
		  'annot_transcript_name', 'n_exons', 'length',\
		  'gene_novelty', 'transcript_novelty',	'ISM_subtype']
fields.extend(datasets)

df = pd.read_csv(infile, sep='\t')  

df = df[fields]

infile = infile.split('talon_abundance')
print(infile)

df.to_csv(infile, index=False, sep='\t')

