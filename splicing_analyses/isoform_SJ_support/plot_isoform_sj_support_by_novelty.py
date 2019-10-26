import pandas as pd
from collections import defaultdict
import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def get_args():

	parser = argparse.ArgumentParser()
	parser.add_argument('-c', dest='csv', 
		help='csv summary file from gen_isoform_support_table.py')
	parser.add_argument('-s', dest='sample_name', 
		help='name of input sample')

	args = parser.parse_args()
	return args

def plot_support(df, args):

	# color palette
	color_dict = defaultdict()
	green = "#009E73"
	orange = "#D55E00"
	gold = "#E69F00"
	blue = "#0072B2"
	pink = "#CC79A7"
	black = "#000000"

	color_dict['Known'] = (green,0)
	color_dict['ISM'] = (blue,1)
	color_dict['NIC'] = (orange,2)
	color_dict['NNC'] = (gold,3)
	color_dict['Intergenic'] = (pink,4)
	color_dict['Antisense'] = (black,5)

	colors = [i[1][0] for i in sorted([(k,color_dict[k]) for k in df['Isoform Novelty'].tolist()], key=lambda x: x[1][1])]
	order = [i[0] for i in sorted([(k,color_dict[k]) for k in df['Isoform Novelty'].tolist()], key=lambda x: x[1][1])]

	# if 'Intergenic' not in df["Isoform Novelty"].tolist():
	# 	colors = [green, blue, orange, gold, black]
	# 	order = ['Known', 'ISM', 'NIC', 'NNC', "Antisense"]
	# else:
	# 	colors = [green, blue, orange, gold, pink, black]
	# 	order = ['Known', 'ISM', 'NIC', 'NNC', 'Intergenic', "Antisense"]

	# plotting
	plt.figure(figsize=(8.5,8.5))
	sns.set(font_scale=1.5, style="whitegrid")

	top_plot = sns.barplot(x='Isoform Novelty', y='total_percent', data=df, 
		color='white', order=order, edgecolor='black')
	bottom_plot = sns.barplot(x='Isoform Novelty', y='percent_illumina', 
		data=df, palette=colors, saturation=1, order=order,
		edgecolor='black')

	topbar = plt.Rectangle((0,0),1,1,fc='white', edgecolor='black')
	bottombar = plt.Rectangle((0,0),1,1,fc='#0000A3',  edgecolor='black')

	plt.title('{} Splice Junction Illumina Support by Isoform Novelty'.format(args.sample_name))
	bottom_plot.set_ylabel("Percent of Isoforms with 100% SJ Support")
	for ntype, p in zip(order, bottom_plot.patches):
		height = p.get_height()
		bottom_plot.text(p.get_x()+p.get_width()/2.,
				height + .3,
				'n = {}'.format(df.loc[df['Isoform Novelty']==ntype]['total'].values[0]),
				ha="center")
	for ntype, p in zip(order, bottom_plot.patches):
		val = df.loc[df['Isoform Novelty'] == ntype]['percent_illumina'].values[0]
		height = val
		height = (height/2)
		bottom_plot.text(p.get_x()+p.get_width()/2.,
				height,
				'{:1.2f}%'.format(val),
				ha="center",
				c='#FFFFFF')

	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_sj_novelty_illumina_isoform_novelty_support.pdf')
	plt.savefig('figures/'+args.sample_name.replace(' ','_')+'_sj_novelty_illumina_isoform_novelty_support.png', dpi = 600)


def main():

	args = get_args()

	df = pd.read_csv(args.csv)
	df.rename({'novelty': 'Isoform Novelty'}, axis=1, inplace=True)
	df.set_index('Isoform Novelty', inplace=True)
	df.rename({'intergenic': 'Intergenic', 'antisense': 'Antisense'}, inplace=True)
	df.reset_index(inplace=True)
	df['total_percent'] = 100


	# print(df)

	plot_support(df, args)

if __name__ == '__main__': main()
