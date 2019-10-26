import pandas as pd 
import argparse
import numpy as np


def get_args():

	desc = 'Determine novelty of SJs in an input file as compared to a ref. '+\
		   'Adds a novelty column to -sj input file.'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-sj', dest='sj_file', 
		help = 'Splice junction file to label SJ novelty')
	parser.add_argument('-ref_sj', dest='ref_sj_file', 
		help = 'Reference splice junction file')
	# parser.add_argument('-sample', dest='sample_name', 
	# 	help = 'Sample name ie GM12878')

	args = parser.parse_args()
	return args

def read_sj_file(infile):
  
	df = pd.read_csv(infile, sep='\t',
		 names=['chrom', 'start', 'stop', 'strand'],
		 usecols=[0,1,2,3])
	df.drop_duplicates(inplace=True)
	return df

def read_extra_info(infile):

	df = pd.read_csv(infile, sep='\t',  
		 names=['chrom', 'start', 'stop', 'strand',
				'motif', 'ann_stat', 'uniq_reads',
				'multi_reads', 'overhang'])
	return df

def find_known_sjs(df, ref_df):

	temp = df.merge(ref_df, how='left', 
		on=['chrom', 'start', 'stop', 'strand'])

	return temp

def get_donors_acceptors(df):

	df['donor'] = df.apply(lambda x: x.start if x.strand==1 else x.stop, axis=1)
	df['acceptor'] = df.apply(lambda x: x.stop if x.strand==1 else x.start, axis=1)

	return df

def main():

	args = get_args()
	sj_file = args.sj_file
	ref_sj_file = args.ref_sj_file

	# get dataframes with whole splice junctions
	sjs = read_sj_file(sj_file)
	ref = read_sj_file(ref_sj_file)
	ref['novelty'] = 'Known'

	# get extra information from the rest of the input file
	sjs_info = read_extra_info(sj_file)

	# label known sjs as such
	sjs = find_known_sjs(sjs, ref)

	# get known sjs
	known_sjs = sjs.copy()
	known_sjs.dropna(inplace=True)
	known_sjs['NNC_type'] = np.NaN

	# subset sjs to get just those that didn't intersect with Gencode
	# these are our candidate NNC or NIC junctions
	novel_sjs = sjs.copy()
	novel_sjs = novel_sjs[novel_sjs[['chrom', 'novelty']].isnull().any(axis=1)]
	
	# create donor and acceptor lists for both the reference and the candidate
	# NNC/NIC junctions
	novel_sjs = get_donors_acceptors(novel_sjs)
	ref = get_donors_acceptors(ref)
	ref_donors = ref.donor.values.tolist()
	ref_acceptors = ref.acceptor.values.tolist()

	# find NICs
	novel_sjs['novelty'] = novel_sjs.apply(
		lambda x: 'NIC' if x.donor in ref_donors and x.acceptor in ref_acceptors else np.NaN,
		axis=1)

	# the rest are NNCs, but does one donor/acceptor match or do none?
	nnc_sjs = novel_sjs.copy()
	nnc_sjs = nnc_sjs[nnc_sjs.isnull().any(axis=1)]
	nnc_sjs['novelty'] = 'NNC'
	nnc_sjs['NNC_type'] = nnc_sjs.apply(
		lambda x: 'NNC_0' if x.donor not in ref_donors and x.acceptor not in ref_acceptors else 'NNC_1',
		axis=1)

	# get the NIC df ready to rumble
	novel_sjs.dropna(inplace=True) # remove NNCs
	novel_sjs['NNC_type'] = np.NaN # set last column to NaN

	# remove acceptor/donor cols
	novel_sjs.drop(['donor', 'acceptor'], axis=1, inplace=True) 
	nnc_sjs.drop(['donor', 'acceptor'], axis=1, inplace=True)
   
	# concatenate known, NIC, and NNC dataframes
	temp = pd.concat([known_sjs, novel_sjs, nnc_sjs])

	# print(temp.head())
	# print(temp.columns)

	# merge temp with its info
	temp = temp.merge(sjs_info, on=['chrom', 'start', 'stop', 'strand'])

	# print(temp.head())
	# print(temp.columns)

	# reorder columns
	temp = temp[['chrom', 'start', 'stop', 'strand',
			'motif', 'ann_stat', 'uniq_reads',
			'multi_reads', 'overhang', 'novelty', 'NNC_type']]

	# write to output file
	ofile = sj_file.replace('.tab', '_novelty.tab') 
	temp.to_csv(ofile, sep='\t', 
		header=False, index=False, na_rep="NaN")

if __name__ == '__main__': main()