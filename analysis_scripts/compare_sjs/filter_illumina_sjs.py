import pandas as pd
import argparse
import numpy as np

def get_args():

	desc = 'Outputs an SJ file with filtered Illumina SJs given sj files from each rep'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-sj_1', dest='sj_1', 
		help = 'Illumina Rep1 SJ file')
	parser.add_argument('-sj_2', dest='sj_2',
		help = 'Illumina Rep2 SJ file')

	return parser.parse_args()

def read_sj_file(infile):
  
	df = pd.read_csv(infile, sep='\t',
		 names=['chrom', 'start', 'stop', 'strand', 'annotated', 'n_reads'],
		 usecols=[0,1,2,3,5,6])
	# df['rep_num'] = rep_num
	# df['sj_id'] = df.apply(lambda x: 
	# 	'_'.join([str(i) for i in [x.chrom, x.start, x.stop, x.strand]]), axis=1)

	return df

def filter_dfs(df1, df2):


	# get rid of everything with 0 supporting reads
	df1 = df1[df1.n_reads != 0]
	print(df1.n_reads.unique()[:10])
	print(df1.head())
	print(len(df1.index))
	df2 = df2[df2.n_reads != 0]
	print(df2.n_reads.unique()[:10])
	print(df2.head())
	print(len(df2.index))

	# merge dfs, preserving read counts for each sj
	df = df1.merge(df2, how='outer', 
		 on=['chrom', 'start', 'stop', 'strand', 'annotated'],
		 suffixes=('_rep1', '_rep2'))
	df = df.fillna(0)

	# remove novel sjs that aren't seen in both reps
	print(len(df.index))
	df.to_csv('before_removing.csv')
	df['reproducible_novel'] = df.apply(lambda x:
		True if x.annotated == False and x.n_reads_rep1 > 0 and x.n_reads_rep2 > 0 else False, axis=1)
	df['reproducible_novel'] = df.apply(lambda x: 
		True if x.annotated == True else x.reproducible_novel, axis=1)
	df = df[df.reproducible_novel == True]
	df.to_csv('after_removing.csv')
	print(len(df.index))

	return df
	# exit()

def main():

	args = get_args()

	df1 = read_sj_file(args.sj_1)
	df2 = read_sj_file(args.sj_2)

	df = filter_dfs(df1, df2)

	df.to_csv(args.sj_1.replace('Rep1_', ''),
		header=False, index=False, sep='\t')

if __name__ == '__main__': main()