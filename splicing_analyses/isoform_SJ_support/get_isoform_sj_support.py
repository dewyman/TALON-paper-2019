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
	if get_field_value('antisense_transcript', fields) == 'TRUE':
		nov_types.append('antisense')
	if get_field_value('intergenic_transcript', fields) == 'TRUE':
		nov_types.append('intergenic')
	if get_field_value('ISM_transcript', fields) == 'TRUE':
		nov_types.append('ISM')
	if get_field_value('NIC_transcript', fields) == 'TRUE':
		nov_types.append('NIC')
	if get_field_value('NNC_transcript', fields) == 'TRUE':
		nov_types.append('NNC')
	if get_field_value('transcript_status', fields) == 'KNOWN':
		nov_types.append('Known')
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

def determine_sj_support(g, df, ref_df, gc_df):

	with open(g, 'r') as infile:

		# loop through gtf
		prev_tid = ""
		prev_exon_end = 0
		n_exons = 0

		for line in infile: 
			line = line.strip().split('\t')

			if line[2] != 'exon': continue
			fields = line[-1]

			tid = get_field_value('transcript_id', fields)
			exon_num = int(get_field_value('exon_number', fields))
			strand = line[6]

			# testing 
			# print(get_field_value('exon_id', fields))

			if tid != prev_tid:

				# print('transcript: '+tid)

				# check if last transcript was monoexonic
				if n_exons == 1: 
					df = df[df.transcript_id != prev_tid]

				# reset number of exons
				n_exons = 1 

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

				# incrememnt number of exons in this transcript
				n_exons += 1

				# get the current junction and check if it's in our list of sjs
				sj = get_sj(line, prev_exon_end)
				# print(sj)

				if strand == "+":
					prev_exon_end = line[4]
				else:
					prev_exon_end = line[3] 

				if sj not in ref_df.id.values.tolist():
					df.loc[df.transcript_id == tid, 'illumina_sj_support'] = False
				if sj not in gc_df.id.values.tolist():
					# print('culprit exon not found in gc~')
					df.loc[df.transcript_id == tid, 'gencode_sj_support'] = False
			# print()
		# clean up if the last entry is monoexonic
		if n_exons == 1: 
				df = df[df.transcript_id != prev_tid]
	
		return df

def main():

	# time execution
	start_time = time.time()

	args = get_args()

	print('Loading illumina sj file....')
	ill_df = read_sj_file(args.ref_sj)
	print('Loading gencode sj file....')
	gc_df = read_sj_file(args.ref_sj_2)

	# # testing
	gc_df.to_csv('gencode_sjs.csv', index=False)

	df = get_transcript_df(args.gtf)
	# print(df.head())
	# print(df.loc[df.transcript_id == 'ENST00000507486.1'])
	# print(df.novelty.unique())

	print('Determining sj support for each transcript....')
	df = determine_sj_support(args.gtf, df, ill_df, gc_df)

	df.to_csv('{}_isoform_sj_support.csv'.format(args.sample_name),
			  index=False)

	# print end time
	print() 
	print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__': main()