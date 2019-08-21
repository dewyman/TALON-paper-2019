import pandas as pd 
import argparse
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import math

def get_args():
	desc = 'Given sj files, see which splice junctions are shared/unique between datasets'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-sj_1', dest='sj_1', 
		help = '1st splice junction file')
	parser.add_argument('-sj_1_name', dest='sj_1_name', 
		help = '1st splice junction file sample name ie "Gencode"')
	parser.add_argument('-sj_2', dest='sj_2', 
		help = '2nd splice junction file')
	parser.add_argument('-sj_2_name', dest='sj_2_name', 
		help = '2nd splice junction file sample name ie "Gencode"')
	parser.add_argument('-sj_3', dest='sj_3',
		help = '3rd splice junction file')
	parser.add_argument('-sj_3_name', dest='sj_3_name', 
		help = '3rd splice junction file sample name ie "Gencode"')
	parser.add_argument('-sample', dest='sample_name', 
		help = 'Sample name ie "PacBio GM12878"')
	parser.add_argument('--log', dest='log_sizes', default=False,
		action='store_true', help = 'Log the sizes of the circles')

	args = parser.parse_args()
	return args

def read_sj_file(infile, dtype):

	df = pd.read_csv(infile, sep='\t', 
		names=['chrom', 'start', 'stop'], usecols=[0,1,2])
	df.drop_duplicates(inplace=True)
	return df

def find_intersect_counts(dfa, dfb, dfc, args):
	
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
	labels = (args.sj_1_name, args.sj_2_name, args.sj_3_name)

	return counts, labels

def main():

	args = get_args()

	# read in each of the sj dfs
	pb_df = read_sj_file(args.sj_1, args.sj_1_name)
	print(len(pb_df.index))
	ont_df = read_sj_file(args.sj_2, args.sj_2_name)
	print(len(ont_df.index))
	ill_df = read_sj_file(args.sj_3, args.sj_3_name)
	print(len(ill_df.index))

	# get each of the intersection counts that we care about
	counts, labels = find_intersect_counts(pb_df, ont_df, ill_df, args)

	# change circle sizes 
	if args.log_sizes:
		intersection_labels = tuple([str(i) for i in counts])
		counts = tuple([math.log2(i) for i in counts])

	print(counts)
	print(labels)

	# plot the venn diagram
	plt.figure(figsize=(8.5,8.5))
	v = venn3(subsets=counts, set_labels=('A', 'B', 'C'))

	# messing with label text
	v.get_label_by_id('A').set_text(args.sj_1_name)
	v.get_label_by_id('B').set_text(args.sj_2_name)
	v.get_label_by_id('C').set_text(args.sj_3_name)
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

	if args.log_sizes:
		v.get_label_by_id('100').set_text(intersection_labels[0])
		v.get_label_by_id('010').set_text(intersection_labels[1])
		v.get_label_by_id('001').set_text(intersection_labels[3])
		v.get_label_by_id('110').set_text(intersection_labels[2])
		v.get_label_by_id('101').set_text(intersection_labels[4])
		v.get_label_by_id('011').set_text(intersection_labels[5])
		v.get_label_by_id('111').set_text(intersection_labels[6])

	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_venn.pdf')
	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_venn.png', dpi = 600)

if __name__ == '__main__':
	main()
