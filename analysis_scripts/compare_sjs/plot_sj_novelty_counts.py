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
	parser.add_argument('--extra_support', dest='extra_file', default=None, 
		help = 'Splice junction file to compare to for extra support')
	parser.add_argument('--support_name', dest='sup_name', default=None,
		help = 'Type of extra support, ie. Illumina')
	parser.add_argument('--nnc_subtypes', dest='split_nnc', default=False, 
		action='store_true', help='Use this flag to plot the NNC subtypes')

	args = parser.parse_args()
	return args

def read_sj_file(infile):
  
	df = pd.read_csv(infile, sep='\t',
		 names=['novelty', 'NNC_type'],
		 usecols=[9,10])
	return df

def read_sjs_sj_file(infile):
  
	df = pd.read_csv(infile, sep='\t',
		 names=['chrom', 'start', 'stop', 'strand', 
		   		'novelty', 'NNC_type'],
		 usecols=[0,1,2,3,9,10])
	return df

def intersect_extra_support(df1, df2):

	# add columns for supported by second set of sjs
	temp = df1.merge(df2, on=['chrom', 'start', 'stop', 'strand',
							  'novelty', 'NNC_type'])

	known_count = len(temp.loc[temp.novelty == 'Known'])
	nic_count = len(temp.loc[temp.novelty == 'NIC'])
	nnc_count = len(temp.loc[temp.novelty == 'NNC'])
	
	counts = [known_count, nic_count, nnc_count]

	return counts

def vanilla_plot(counts, args):

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
					 palette=colors, saturation=1, order=order, edgecolor='black')
	plt.ylim(0, 20)
	ax.set_ylabel("log2(SJ count)")
	plt.title(args.sample_name+' Splice Junction Novelty')
	for ntype, p in zip(order, ax.patches):
		height = p.get_height()
		ax.text(p.get_x()+p.get_width()/2.,
				height + .3,
				'{:1.2f}%'.format(nov_counts.loc[nov_counts.Novelty == ntype]['percent'].values[0]),
				ha="center")

	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_sj_novelty.pdf')
	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_sj_novelty.png', dpi = 600)

def extra_support_plot(counts, sup_counts, args):

	# get no. unsupported sjs 
	unsup_counts = [counts[i] - sup_counts[i] for i in range(len(counts))]

	# collapse df to one we can use with seaborn
	nov_counts = pd.DataFrame(data=counts, columns=['Count'], 
		index=['Known', 'NIC', 'NNC'])
	nov_counts.reset_index(inplace=True)
	nov_counts.rename({'index': 'Novelty'}, inplace=True, axis=1)
	nov_counts['Supported Count'] = sup_counts
	nov_counts['Unsupported Count'] = unsup_counts

	# add aesthetics-related columns
	nov_counts['log2(count)'] = nov_counts.apply(
		lambda x: math.log2(x.Count), axis=1)
	nov_counts['percent'] = nov_counts.apply(
		lambda x: 100*(float(x['Supported Count'])/x.Count), axis=1)
	nov_counts['log2(supported count)'] = nov_counts.apply(
		lambda x: (x.percent/100)*x['log2(count)'], axis=1)
	nov_counts['total_percent'] = 100

	# color palette
	green = "#009E73"
	orange = "#D55E00"
	gold = "#E69F00"
	colors = [green, orange, gold]

	# plotting
	plt.figure(figsize=(8.5,8.5))
	sns.set(font_scale=1.5, style="whitegrid")
	order = ['Known', 'NIC', 'NNC']

	top_plot = sns.barplot(x='Novelty', y='total_percent', data=nov_counts, 
		color='white', order=order, edgecolor='black')
	bottom_plot = sns.barplot(x='Novelty', y='percent', 
		data=nov_counts, palette=colors, saturation=1, order=order,
		edgecolor='black')

	topbar = plt.Rectangle((0,0),1,1,fc='white', edgecolor='black')
	bottombar = plt.Rectangle((0,0),1,1,fc='#0000A3',  edgecolor='black')

	# plt.ylim(0, 20)
	plt.title('{} Splice Junction Novelty with {} Support'.format(args.sample_name,
		args.sup_name))
	bottom_plot.set_ylabel("Percent of Total SJs")
	for ntype, p in zip(order, bottom_plot.patches):
		height = p.get_height()
		bottom_plot.text(p.get_x()+p.get_width()/2.,
				height + .3,
				'n = {}'.format(nov_counts.loc[nov_counts.Novelty==ntype]['Count'].values[0]),
				ha="center")
	for ntype, p in zip(order, bottom_plot.patches):
		val = nov_counts.loc[nov_counts.Novelty == ntype]['percent'].values[0]
		height = val
		height = (height/2)
		bottom_plot.text(p.get_x()+p.get_width()/2.,
				height,
				'{:1.2f}%'.format(val),
				ha="center")

	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_sj_novelty_{}_support.pdf'.format(args.sup_name))
	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_sj_novelty_{}_support.png'.format(args.sup_name), dpi = 600)

def main():

	args = get_args()

	df = read_sj_file(args.sj_file)

	# extra support 
	if args.extra_file:
		main_df = read_sjs_sj_file(args.sj_file)
		extra_df = read_sjs_sj_file(args.extra_file)

		# get counts supported by extra sj file
		sup_counts = intersect_extra_support(main_df, extra_df)

	# get known, NIC, and NNC counts
	known_count = len(df.loc[df.novelty == 'Known'].index)
	nic_count = len(df.loc[df.novelty == 'NIC'].index)
	nnc_count = len(df.loc[df.novelty == 'NNC'].index)
	counts = [known_count, nic_count, nnc_count]

	# plot whichever thing we want
	if args.extra_file: extra_support_plot(counts, sup_counts, args)
	else: vanilla_plot(counts, args)

if __name__ == '__main__': main()