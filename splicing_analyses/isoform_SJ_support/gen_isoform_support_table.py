import pandas as pd 
import argparse

def get_args():

	desc = 'Get a summary table of isoform SJ support by transcript novelty type'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-csv', dest='support_csv',
		help = 'CSV of SJ support output from get_isoform_support.py')
	parser.add_argument('-sample', dest='sample_name', 
		help = 'Sample name ie GM12878')

	args = parser.parse_args()
	return args



def read_support_file(csvfile):

	df = pd.read_csv(csvfile)
	df.drop('transcript_id', inplace=True, axis=1)
	return df

def gen_summary_table(df):

	# groupby novelty and whichever support and count
	ill_sup = df.groupby(['novelty','illumina_sj_support']).count()
	ill_sup.reset_index(inplace=True)
	ill_sup['total'] = ill_sup.apply(lambda x:
			len(df.loc[df.novelty == x.novelty].index), axis=1)
	ill_sup.rename({'illumina_sj_support': 'support', 
	               'gencode_sj_support': 'illumina_count'},
	                inplace=True, axis=1)
	ill_sup = ill_sup[ill_sup.support != False]
	ill_sup['percent_illumina'] = ill_sup.apply(lambda x: 
	        (x.illumina_count/x.total)*100, axis=1)
	ill_sup.drop('support', inplace=True, axis=1)

	gc_sup = df.groupby(['novelty','gencode_sj_support']).count()
	gc_sup.reset_index(inplace=True)
	gc_sup['total'] = gc_sup.apply(lambda x:
		len(df.loc[df.novelty == x.novelty].index), axis=1)
	gc_sup.rename({'gencode_sj_support': 'support', 
	           'illumina_sj_support': 'gencode_count'},
	            inplace=True, axis=1)
	gc_sup = gc_sup[gc_sup.support != False]
	gc_sup['percent_gencode'] = gc_sup.apply(lambda x: 
	    (x.gencode_count/x.total)*100, axis=1)
	gc_sup.drop('support', inplace=True, axis=1)

	# furnish out the final table
	ill_sup = ill_sup.merge(gc_sup, how='outer', on=['novelty','total'])
	ill_sup = ill_sup[['novelty', 'total', 'illumina_count', 
					   'percent_illumina', 'gencode_count', 
					   'percent_gencode']]
	ill_sup.fillna(0, inplace=True)

	return ill_sup

def main():

	args = get_args()

	df = read_support_file(args.support_csv)

	df = gen_summary_table(df)

	df.to_csv('{}_isoform_sj_support_summary.csv'.format(args.sample_name),
	          index=False)


if __name__ == '__main__': main()