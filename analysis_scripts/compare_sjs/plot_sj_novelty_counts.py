import pandas as pd 
import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import math 

def get_args():

	desc = 'Determine novelty of SJs in an input file as compared to a ref. '+\
		   'Adds a novelty column to -sj input file.'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-sj', dest='sj_file', 
		help = 'Splice junction file to label SJ novelty')
	parser.add_argument('-sample', dest='sample_name', 
		help = 'Sample name ie GM12878')

	args = parser.parse_args()
	return args

def read_sj_file(infile):
  
	df = pd.read_csv(infile, sep='\t',
		 names=['novelty', 'NNC_type'],
		 usecols=[9,10])
	return df

def main():

	args = get_args()

	df = read_sj_file(args.sj_file)

	# get known, NIC, and NNC counts
	known_count = len(df.loc[df.novelty == 'Known'].index)
	nic_count = len(df.loc[df.novelty == 'NIC'].index)
	nnc_count = len(df.loc[df.novelty == 'NNC'].index)
	counts = [known_count, nic_count, nnc_count]

	# collapse df to one we can use with seaborn
	nov_counts = pd.DataFrame(data=counts, columns=['Count'], 
		index=['Known', 'NIC', 'NNC'])
	nov_counts.reset_index(inplace=True)
	nov_counts.rename({'index': 'Novelty'}, inplace=True, axis=1)

	# add aesthetics-related columns
	nov_counts['log2(count)'] = nov_counts.apply(
		lambda x: math.log2(x.Count), axis=1)
	nov_counts['percent'] = nov_counts.apply(
		lambda x: 100*(float(x.Count)/sum(counts)), axis=1)
	# add illumina support/unsupport columns here to use for stacked bar later!

	# color palette
	green = "#009E73"
	orange = "#D55E00"
	gold = "#E69F00"
	colors = [green, orange, gold]

	# plotting
	plt.figure(figsize=(8.5,8.5))
	sns.set(font_scale=1.5, style="whitegrid")
	order = ['Known', 'NIC', 'NNC']
	ax = sns.barplot(x='Novelty', y='log2(count)', data=nov_counts,
					 palette=colors, saturation=1, order=order)
	plt.ylim(0, 20)
	plt.title(args.sample_name+' Splice Junction Novelty')
	for ntype, p in zip(order, ax.patches):
		height = p.get_height()
		ax.text(p.get_x()+p.get_width()/2.,
				height + .3,
				'{:1.2f}%'.format(nov_counts.loc[nov_counts.Novelty == ntype]['percent'].values[0]),
				ha="center")

	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_sj_novelty.pdf')
	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_sj_novelty.png', dpi = 600)

if __name__ == '__main__': main()