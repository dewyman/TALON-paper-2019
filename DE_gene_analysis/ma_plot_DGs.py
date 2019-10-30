import seaborn as sns
import pandas as pd
import argparse

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

	df = pd.read_csv(f, sep='\t',
		names=['gene', 'log2FC', 'log2CPM', 'p'],
		usecols=[0,1,2,4],
		header=0)
	return df

def plot_pairwise_MA(dfa, dfb, namea, nameb):

	# merge on gene name
	df = pd.merge(dfa, dfb, on='gene', how='inner', suffixes=('_a', '_b'), sort=True)

	# assign labels based on whether things are jointly differentially expressed or nah
	df['de_category'] = df.apply(lambda x: 'NS' if x.p_a == 1 and x.p_b == 1 else 
						('DE_both' if x.p_a != 1 and x.p_b != 1 else
						('DE_a' if x.p_a != 1 else 'DE_b')), axis=1)

	# MA plot colored by which technologies identified each gene as DE
	p = sns.scatterplot(data=df, x='log2FC_a', y='log2CPM_a', hue='de_category')
	fig = p.get_figure()
	fig.savefig('test_pb.png')

	p = sns.scatterplot(data=df, x='log2FC_b', y='log2CPM_b', hue='de_category')
	fig = p.get_figure()
	fig.savefig('test_illumina.png')

	exit()

def main():
	args = get_args()

	pb_df = read_dge_file(args.pbfile)
	# print(pb_df.dtypes)
	ont_df = read_dge_file(args.ontfile)
	# print(ont_df.dtypes)
	ill_df = read_dge_file(args.illfile)
	# print(ill_df.dtypes)

	plot_pairwise_MA(pb_df, ill_df, 'PacBio', 'Illumina')

if __name__ == '__main__': main()
