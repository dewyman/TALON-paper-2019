import pandas as pd 
import argparse

parser = argparse.ArgumentParser(description=\
		  'From an exisiting abundance file, cut out\
		   input dataset fields')
parser.add_argument('--infile', help='input abundance file')
parser.add_argument('--outfile', help='output abundance file', default=False)
parser.add_argument('--d', help='list of datasets to include, one on each line')
args = parser.parse_args()

infile = args.infile
if args.outfile:
	ofile = args.outfile
else:
	ofile = infile

dfile = open(args.d, 'r')
datasets = []
for line in dfile:
	line = line.replace('\n', '')
	datasets.append(line)

fields = ['gene_ID','transcript_ID','annot_gene_id',\
		  'annot_transcript_id', 'annot_gene_name', \
		  'annot_transcript_name', 'n_exons', 'length',\
		  'gene_novelty', 'transcript_novelty',	'ISM_subtype']
fields.extend(datasets)

df = pd.read_csv(infile, sep='\t')  

df = df[fields]

df.to_csv(ofile, index=False, sep='\t')

