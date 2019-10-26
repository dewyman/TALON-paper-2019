import os
import pandas as pd
import glob
import seaborn as sns
import argparse
import time


def get_args():

	desc = 'Plots the support for transcript isoforms based on Illumina data'
	parser = argparse.ArgumentParser(description=desc)

	parser.add_argument('-gtf', dest='gtf', 
		help = 'GTF from TALON')
	parser.add_argument('-ref_sj_1', dest='ref_sj', 
		help = 'Reference SJ file (Illumina)')
	parser.add_argument('-ref_sj_2', dest='ref_sj_2')
	parser.add_argument('-sample', dest='sample_name', 
		help = 'Sample name ie GM12878')

	return parser.parse_args()

def read_sj_file(infile):
  
	df = pd.read_csv(infile, sep='\t',
		 names=['chrom', 'start', 'stop', 'strand'],
		 usecols=[0,1,2,3])
	df['id'] = df.apply(
		lambda x: '_'.join([str(i) for i in [x.chrom, x.start, x.stop, x.strand]]),
		axis = 1)
	df.drop(['chrom', 'start', 'stop', 'strand'], inplace=True, axis=1)

	return df

# get value associated with keyword in the 9th column of gtf
def get_field_value(key, fields):
	if key not in fields:
		return None
	else:
		return fields.split(key+' "')[1].split()[0].replace('";','')

def find_novelty_type(fields):
	nov_types = []
	if 'antisense' in fields: nov_types.append('antisense')
	if 'intergenic' in fields: nov_types.append('intergenic')
	if 'ISM_transcript' in fields: nov_types.append('ISM')
	if 'NIC' in fields: nov_types.append('NIC')
	if 'NNC' in fields: nov_types.append('NNC')
	if 'KNOWN' in fields: nov_types.append('Known')
	return nov_types

def get_transcript_df(g):

	df = pd.DataFrame(columns=['transcript_id', 'novelty'])
	with open(g, 'r') as infile:

		for line in infile:
			line = line.strip().split('\t')
			fields = line[-1]

			if line[2] != 'transcript': continue

			tid = get_field_value('transcript_id', fields)

			# find novelty type(s)
			nov_types = find_novelty_type(fields)
			# if len(nov_types) > 1: 
			# 	print(nov_types)
			# 	print(tid)
			for n in nov_types: 
				df = df.append({'transcript_id': tid, 'novelty': n}, 
					ignore_index=True)

	# add sj support column
	df['illumina_sj_support'] = True
	df['gencode_sj_support'] = True

	return df

def get_sj(entry, prev_exon_end, minIntron=21):

	chromosome = entry[0]
	strand = entry[6]

	if strand == "+": 
		strand = "1"
		intron_start = int(prev_exon_end) + 1
		intron_end = int(entry[3]) - 1

	elif strand == "-":
		strand = "2"
		intron_start = int(entry[4]) + 1  
		intron_end = int(prev_exon_end) - 1	
	
	if abs(intron_end - intron_start + 1) < minIntron:
		return None 
	
	return "_".join([chromosome, str(intron_start), str(intron_end), strand])

def determine_sj_novelty(g, df, ref_df, gc_df):

	with open(g, 'r') as infile:

		# loop through gtf
		prev_tid = ""
		prev_exon_end = 0
		for line in infile: 
			line = line.strip().split('\t')

			if line[2] != 'exon': continue
			fields = line[-1]

			tid = get_field_value('transcript_id', fields)
			exon_num = int(get_field_value('exon_number', fields))
			strand = line[6]

			if tid != prev_tid:

				# new transcript
				if exon_num != 1:
					print("Error: exons are not listed in order")
					exit()
				prev_tid = tid
				if strand == "+":
					prev_exon_end = line[4]
				else:
					prev_exon_end = line[3]
			else: 
				# get the current junction and check if it's in our list of sjs

				sj = get_sj(line, prev_exon_end)

				if strand == "+":
					prev_exon_end = line[4]
				else:
					prev_exon_end = line[3] 

				if sj not in ref_df.id.values.tolist():
					df.loc[df.transcript_id == tid, 'illumina_sj_support'] = False
				if sj not in gc_df.id.values.tolist():
					df.loc[df.transcript_id == tid, 'gencode_sj_support'] = False

	return df

def main():
	
	# time execution
	start_time = time.time()

	args = get_args()

	ill_df = read_sj_file(args.ref_sj)
	gc_df = read_sj_file(args.ref_sj_2)

	df = get_transcript_df(args.gtf)
	# print(df.head())
	# print(df.loc[df.transcript_id == 'ENST00000507486.1'])
	# print(df.novelty.unique())

	df = determine_sj_novelty(args.gtf, df, ill_df, gc_df)

	df.to_csv('isoform_sj_support.csv', index=False)

	# print end time
	print()
	print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__': main()