import seaborn as sns
import pandas as pd
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt

def get_args():
	desc = 'Given sj files, see which splice junctions are shared/unique between datasets'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-pb', dest='pbfile', 
		help = 'PacBio DE genes table')
	parser.add_argument('-ont', dest='ontfile', 
		help = 'ONT DE genes table')
	parser.add_argument('-ill', dest='illfile',
		help = 'Illumina DE genes table')
	parser.add_argument('-sample', dest='sample_name', 
		help = 'Sample name "HepG2 vs. K562"')

	args = parser.parse_args()
	return args

def read_dge_file(f):

	try:
		df = pd.read_csv(f, sep='\t',
			names=['gene', 'log2FC', 'log2CPM', 'p'],
			usecols=[0,1,2,5],
			header=0)
	except:
		df = pd.read_csv(f, sep=',',
			names=['gene', 'log2FC', 'log2CPM', 'p'],
			usecols=[0,1,2,5],
			header=0)
	return df

def plot_pairwise_MA(dfa, dfb, namea, nameb, sample_name):

	# color palette
	color_dict = defaultdict()
	blue = '#999cfc'
	red = '#fd9a9b'
	red3 = "#fc0509"
	green = '#9acb9b'
	yellow = '#fda429'
	magenta = '#df9bde'
	gray ='#999999'
	green2 = "#009E73"
	blue2 = "#56B4E9"
	blue5 = "#87CEEB"
	blue3 = "#0072B2"
	red2 = "#D55E00" 
	pink = "#CC79A7"
	blue4 = "#000080"
	color_dict['DE in {}'.format(namea)] = red3
	color_dict['DE in {}'.format(nameb)] = blue5
	color_dict['DE in both'] = blue4
	color_dict['Not Significant'] = yellow

	# order = []

	# merge on gene name
	df = pd.merge(dfa, dfb, on='gene', how='inner', suffixes=('_a', '_b'), sort=True)

	# assign labels based on whether things are jointly differentially expressed or nah
	df['de_category'] = df.apply(lambda x: 'Not Significant' if x.p_a == 1 and x.p_b == 1 else 
						('DE in both' if x.p_a != 1 and x.p_b != 1 else
						('DE in {}'.format(namea) if x.p_a != 1 else 'DE in {}'.format(nameb))),
						 axis=1)
	df['de_binary'] = df.apply(lambda x: 1 if x.de_category != 'Not Significant' else 0, axis=1)
	# print(df.head())

	# # MA plot colored by which technologies identified each gene as DE
	# plt.figure(figsize=(8.5,8.5))
	# plt.rc('font', size=24)
	# sns.set_style("whitegrid")
	# ax = sns.scatterplot(data=df, x='log2FC_b', y='log2CPM_b', hue='de_category', 
	# 	alpha=0.45, edgecolor=None, size='de_binary', 
	# 	palette=color_dict, sizes={1: 20, 0: 10},
	# 	hue_order=['DE in both', 'DE in {}'.format(namea), 'DE in {}'.format(nameb), 'Not Significant'])
	# plt.xlabel('{} log2-fold change in {}'.format(sample_name, nameb))
	# plt.ylabel('log(Counts per million) in {}'.format(nameb))
	# sns.set(font_scale=1.5, style='whitegrid')
	# plt.xlim(-35, 24)
	# # messing with the legend
	# # removes the title
	# handles, labels = ax.get_legend_handles_labels()
	# ax.legend(loc='lower left', handles=handles[1:-3],
	# 		  labels=labels[1:-3], fancybox=True, framealpha=1)

	# fig = ax.get_figure()
	# fig.savefig('DE_genes_{}_{}_{}_expression.png'.format(sample_name.replace(' ', '_'), namea, nameb))

	# # MA plot colored by which technologies identified each gene as DE
	# plt.figure(figsize=(8.5,8.5))
	# # plt.rc('font', size=24)
	# sns.set_style("whitegrid")
	# ax = sns.scatterplot(data=df, x='log2FC_a', y='log2CPM_a', hue='de_category', 
	# 	alpha=0.45, edgecolor=None, size='de_binary', 
	# 	palette=color_dict, sizes={1: 20, 0: 10},
	# 	hue_order=['DE in both', 'DE in {}'.format(namea), 'DE in {}'.format(nameb), 'Not Significant'])
	# plt.xlabel('{} log2-fold change in {}'.format(sample_name, namea))
	# plt.ylabel('log(Counts per million) in {}'.format(namea))
	# sns.set(font_scale=1.5, style='whitegrid')
	# # plt.xlim(-22, 17)
	# plt.xlim(-35, 24)
	# # messing with the legend
	# # removes the title
	# handles, labels = ax.get_legend_handles_labels()
	# ax.legend(loc='lower left', handles=handles[1:-3],
	# 		  labels=labels[1:-3], fancybox=True, framealpha=1)

	# fig = ax.get_figure()
	# fig.savefig('DE_genes_{}_{}_{}_expression.png'.format(sample_name.replace(' ', '_'), nameb, namea))
	# plot_a_vs_b(df, 'log2FC_a', 'log2CPM_a', namea, nameb, sample_name, color_dict)
	# plot_a_vs_b(df, 'log2FC_b', 'log2CPM_b', nameb, namea, sample_name, color_dict)


	# MA plot colored by which technologies identified each gene as DE
	plt.figure(figsize=(8.5,8.5))
	plt.rc('font', size=24)
	sns.set_style("whitegrid")
	ax = sns.scatterplot(data=df, x='log2FC_a', y='log2CPM_a', hue='de_category', 
		alpha=0.45, edgecolor=None, size='de_binary', 
		palette=color_dict, sizes={1: 20, 0: 10},
		hue_order=['DE in both', 'DE in {}'.format(namea), 'DE in {}'.format(nameb), 'Not Significant'])
	plt.xlabel('{} log2-fold change in {}'.format(sample_name, namea))
	plt.ylabel('log(Counts per million) in {}'.format(namea))
	# sns.set(font_scale=1.5, style='whitegrid')
	# plt.xlim(-22, 17)
	plt.xlim(-35, 25)
	plt.ylim(-1,15)
	# messing with the legend
	# removes the title
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(loc='lower left', handles=handles[1:-3],
			  labels=labels[1:-3], framealpha=1,
			  handletextpad=0.05,
			  fontsize='small',
			  markerscale=1.5,
			  labelspacing=0.1)

	fig = ax.get_figure()
	fig.savefig('DE_genes_{}_{}_{}_expression.png'.format(sample_name.replace(' ', '_'), nameb, namea))

		# MA plot colored by which technologies identified each gene as DE
	plt.figure(figsize=(8.5,8.5))
	plt.rc('font', size=24)
	sns.set_style("whitegrid")
	ax = sns.scatterplot(data=df, x='log2FC_b', y='log2CPM_b', hue='de_category', 
		alpha=0.45, edgecolor=None, size='de_binary', 
		palette=color_dict, sizes={1: 20, 0: 10},
		hue_order=['DE in both', 'DE in {}'.format(namea), 'DE in {}'.format(nameb), 'Not Significant'])
	plt.xlabel('{} log2-fold change in {}'.format(sample_name, nameb))
	plt.ylabel('log(Counts per million) in {}'.format(nameb))
	# sns.set(font_scale=1.5, style='whitegrid')
	# plt.xlim(-22, 17)
	plt.xlim(-35, 25)
	plt.ylim(-1,15)
	# messing with the legend
	# removes the title
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(loc='lower left', handles=handles[1:-3],
			  labels=labels[1:-3], framealpha=1,
			  handletextpad=0.05,
			  fontsize='small',
			  markerscale=1.5,
			  labelspacing=0.1)

	fig = ax.get_figure()
	fig.savefig('DE_genes_{}_{}_{}_expression.png'.format(sample_name.replace(' ', '_'), namea, nameb))

# def plot_a_vs_b(df, field_x, field_y, namea, nameb, sample_name, color_dict):	

	# # MA plot colored by which technologies identified each gene as DE
	# plt.figure(figsize=(8.5,8.5))
	# plt.rc('font', size=24)
	# sns.set_style("whitegrid")
	# ax = sns.scatterplot(data=df, x=field_x, y=field_y, hue='de_category', 
	# 	alpha=0.45, edgecolor=None, size='de_binary', 
	# 	palette=color_dict, sizes={1: 20, 0: 10},
	# 	hue_order=['DE in both', 'DE in {}'.format(namea), 'DE in {}'.format(nameb), 'Not Significant'])
	# plt.xlabel('{} log2-fold change in {}'.format(sample_name, namea))
	# plt.ylabel('log(Counts per million) in {}'.format(namea))
	# # sns.set(font_scale=1.5, style='whitegrid')
	# # plt.xlim(-22, 17)
	# plt.xlim(-35, 25)
	# plt.ylim(-1,15)
	# # messing with the legend
	# # removes the title
	# handles, labels = ax.get_legend_handles_labels()
	# ax.legend(loc='lower left', handles=handles[1:-3],
	# 		  labels=labels[1:-3], framealpha=1,
	# 		  handletextpad=0.05,
	# 		  fontsize='small',
	# 		  markerscale=1.5,
	# 		  labelspacing=0.1)

	# fig = ax.get_figure()
	# fig.savefig('DE_genes_{}_{}_{}_expression.png'.format(sample_name.replace(' ', '_'), nameb, namea))

def main():
	args = get_args()

	pb_df = read_dge_file(args.pbfile)
	# print(pb_df.dtypes)
	ont_df = read_dge_file(args.ontfile)
	# print(ont_df.dtypes)
	ill_df = read_dge_file(args.illfile)
	# print(ill_df.dtypes)

	plot_pairwise_MA(pb_df, ill_df, 'PacBio', 'Illumina', args.sample_name)
	plot_pairwise_MA(ont_df, ill_df, 'ONT', 'Illumina', args.sample_name)

if __name__ == '__main__': main()
