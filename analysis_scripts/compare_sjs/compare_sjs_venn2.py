import pandas as pd 
import argparse
from matplotlib_venn import venn2
import matplotlib.pyplot as plt
import math

def get_args():
	desc = 'Given sj files, see which splice junctions are shared/unique between datasets. Also save the unsupported junctions'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-sj_1', dest='sj_1', 
		help = '1st splice junction file')
	parser.add_argument('-sj_1_name', dest='sj_1_name', 
		help = '1st splice junction file sample name ie "Gencode"')
	parser.add_argument('-sj_2', dest='sj_2', 
		help = '2nd splice junction file')
	parser.add_argument('-sj_2_name', dest='sj_2_name', 
		help = '2nd splice junction file sample name ie "Gencode"')
	parser.add_argument('-sample', dest='sample_name', 
		help = 'Sample name ie "PacBio GM12878"')
	parser.add_argument('--log', dest='log_sizes', default=False,
		action='store_true', help = 'Log the sizes of the circles')

	args = parser.parse_args()
	return args

def find_intersect_counts(dfa, dfb):

	temp = pd.merge(dfa, dfb, how='inner', 
					on=['chrom', 'start', 'stop', 'strand'])
	count_ab = len(temp.index)

	count_a = len(dfa.index) - count_ab
	count_b = len(dfb.index) - count_ab

	counts = (count_a, count_b, count_ab)

	# get the unsupported long-read dfa stuff
	# new = dfa[(~dfa.start.isin(temp.start))&(~dfa.stop.isin(temp.stop))&(~dfa.chrom.isin(temp.chrom))&(~dfa.strand.isin(temp.strand))]
	new = dfa[~dfa.isin(temp)].dropna(how = 'all')
	print('number of unsupported long read sjs : '+str(len(new.index)))

	return (counts,new) 

def read_sj_file(infile):

	df = pd.read_csv(infile, sep='\t', 
		names=['chrom', 'start', 'stop', 'strand'], usecols=[0,1,2,3])
	df.drop_duplicates(inplace=True)
	return df

def main():

	args = get_args()

	dfa = read_sj_file(args.sj_1)
	dfb = read_sj_file(args.sj_2)

	(counts, temp) = find_intersect_counts(dfa, dfb)

	temp.to_csv('{}_unsupported_sjs.tab'.format(args.sample_name).replace(' ', '_'),
		index=False, sep='\t', header=False)

	plt.figure(figsize=(8.5,8.5))
	v = venn2(subsets=counts, set_labels=('A', 'B'))
	plt.title('{} Splice Junction Support'.format(args.sample_name), 
			  fontsize='xx-large')

	# messing with label text
	v.get_label_by_id('A').set_text(args.sj_1_name)
	v.get_label_by_id('B').set_text(args.sj_2_name)
	v.get_label_by_id('A').set_fontsize('x-large')
	v.get_label_by_id('B').set_fontsize('x-large')

	# messing with numerical text
	v.get_label_by_id('10').set_fontsize('x-large')
	v.get_label_by_id('01').set_fontsize('x-large')
	v.get_label_by_id('11').set_fontsize('x-large')

	# save figures
	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_gc_venn2.pdf')
	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_gc_venn2.png',
	            dpi = 600)


if __name__ == '__main__': main()