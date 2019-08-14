import pandas as pd 
import argparse
from matplotlib_venn import venn3
import matplotlib.pyplot as plt

def get_args():
	desc = 'Given sj files, see which splice junctions are shared/unique between datasets'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-pb', dest='pb_sj_file', 
		help = 'Splice junction file from PacBio')
	parser.add_argument('-ont', dest='ont_sj_file', 
		help = 'Splice junction file from ONT')
	parser.add_argument('-illumina', dest='ill_sj_file',
		help = 'Splice junction file from illumina')
	parser.add_argument('-sample', dest='sample_name', 
		help = 'Sample name ie GM12878')

	args = parser.parse_args()
	return args

def read_sj_file(infile, dtype):

	df = pd.read_csv(infile, sep='\t', 
		names=['chrom', 'start', 'stop'], usecols=[0,1,2])
	# df['dtype'] = dtype
	return df

def find_intersect_counts(dfa, dfb, dfc):
	
	# intersection of all (a,b,c) 
	temp = pd.merge(dfa, dfb, how='inner', on=['chrom', 'start', 'stop'])
	temp = pd.merge(temp, dfc, how='inner', on=['chrom', 'start', 'stop'])
	count_abc = len(temp.index)

	# intersection of (a,b)
	temp = pd.merge(dfa, dfb, how='inner', on=['chrom', 'start', 'stop'])
	count_ab = len(temp.index) - count_abc

	# intersection of (a,c)
	temp = pd.merge(dfa, dfc, how='inner', on=['chrom', 'start', 'stop'])
	count_ac = len(temp.index) - count_abc

	# intersection of (b,c)
	temp = pd.merge(dfb, dfc, how='inner', on=['chrom', 'start', 'stop'])
	count_bc = len(temp.index) - count_abc

	# just a
	count_a = len(dfa.index) - count_ab - count_ac - count_abc

	# just b
	count_b = len(dfb.index) - count_ab - count_bc - count_abc

	# just c
	count_c = len(dfc.index) - count_ac - count_bc - count_abc

	counts = (count_a, count_b, count_ab, count_c, count_ac, count_bc, count_abc) 
	labels = ('PacBio', 'ONT', 'Illumina')

	return counts, labels

def main():

	args = get_args()

	# read in each of the sj dfs
	pb_df = read_sj_file(args.pb_sj_file, 'PacBio')
	print(len(pb_df.index))
	ont_df = read_sj_file(args.ont_sj_file, 'ONT')
	print(len(ont_df.index))
	ill_df = read_sj_file(args.ill_sj_file, 'Illumina')
	print(len(ill_df.index))

	# get each of the intersection counts that we care about
	counts, labels = find_intersect_counts(pb_df, ont_df, ill_df)

	print(counts)
	print(labels)

	# plot the venn diagram
	plt.figure(figsize=(8.5,8.5))
	v = venn3(subsets=counts, set_labels=('A', 'B', 'C'))

	# messing with label text
	v.get_label_by_id('A').set_text('PacBio')
	v.get_label_by_id('B').set_text('ONT')
	v.get_label_by_id('C').set_text('Illumina')
	v.get_label_by_id('A').set_fontsize('x-large')
	v.get_label_by_id('B').set_fontsize('x-large')
	v.get_label_by_id('C').set_fontsize('x-large')
	plt.title('{} Splice Junction Support'.format(args.sample_name), fontsize='xx-large')

	# messing with numerical text
	v.get_label_by_id('100').set_fontsize('x-large')
	v.get_label_by_id('010').set_fontsize('x-large')
	v.get_label_by_id('001').set_fontsize('x-large')
	v.get_label_by_id('110').set_fontsize('x-large')
	v.get_label_by_id('101').set_fontsize('x-large')
	v.get_label_by_id('011').set_fontsize('x-large')
	v.get_label_by_id('111').set_fontsize('x-large')

	plt.savefig(args.sample_name+'_venn.pdf')

if __name__ == '__main__':
	main()