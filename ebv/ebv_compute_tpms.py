# get human_ebv_tpms.py
import pandas as pd
import argparse
import os
import math
import datetime
import subprocess

# get basename from a file and path string
def get_basename(filepath):
    import os
    return os.path.basename(os.path.splitext(filepath)[0])    

# get and format output directory
def format_odir(odir):
    import os
    cwd = os.getcwd()

    if odir != '':
        # if first character is not /, use cwd to make this an absolute path
        if odir[0] != '/' and odir[0] != '~':
            odir = cwd+odir
        if odir[-1] != '/':
            odir += '/'
    return odir

# make a dated output directory for the files used for the tracks
def make_dated_folder(odir, bname):
    date = datetime.datetime.now()
    date = date.strftime('%y%m%d')
    odir = odir+date+'_'+bname+'_figures/'

    if not os.path.isdir(odir):
        print('Making output directory '+odir)
        os.makedirs(odir)
    return odir

# get value associated with keyword in the 9th column of gtf
def get_field_value(key, fields):
    if key not in fields:
        return None
    else:
        return fields.split(key+' "')[1].split()[0].replace('";','')

# calculate tpm for a column in the abundance table
def get_tpm(df, col):
	new_col = 'TPM_'+col
	total_reads = df[d].sum()
	df[new_col] = df.apply(lambda x: float(x[d]*1000000)/total_reads, axis=1)

	return new_col, df

# calculate tpm for a column in the abundance table
def get_log_tpm(df, col, gene):
	tpm_col = 'TPM_'+col
	if not gene:
		new_col = 'log_'+tpm_col
	else:
		new_col = 'gene_log_'+TPM_col
	df[new_col] = df.apply(lambda x: math.log2(x[tpm_col]+1), axis=1)

	return new_col, df

# get gtf file name 
parser = argparse.ArgumentParser(description='removes EBV transcripts from GTF file')
parser.add_argument('--c', help='csv config file with the following fields'\
								'1: Full GTF'\
								'2: Full, filtered abundance table'\
                                '3: Full, unfiltered abundance table'\
								'4: EBV, filtered abundance table'\
                                '5: EBV, unfiltered abundance table')
args = parser.parse_args()
config = args.c 
with open(config, 'r') as c:
	for line in c:
		1
line = line.replace('\n', '')
line = line.split(',')
full_gtf = line[0]
full_ab = line[1]
full_unf_ab = line[2]
ebv_ab = line[3]
ebv_unf_ab = line[4]

print(full_gtf)
print(full_ab)
print(full_unf_ab)
print(ebv_ab)
print(ebv_unf_ab)

# get all human transcript ids
infile = open(full_gtf, 'r')
human_tids = []
ebv_tids = []
for i, line in enumerate(infile): 
    line = line.replace('\n', '')
    temp = line.split('\t')
    fields = temp[-1]

    if temp[0] != 'chrEBV' and temp[2] == 'transcript':
        human_tids.append(get_field_value('talon_transcript', fields))
    elif temp[0] == 'chrEBV' and temp[2] == 'transcript':
        ebv_tids.append(get_field_value('talon_transcript', fields))
        
full_df = pd.read_csv(full_ab, sep='\t')
ebv_df = pd.read_csv(ebv_ab, sep='\t')

# reformat human table
# dc_datasets = ['D4', 'D5', 'D10', 'D11']
datasets = ['PacBio_GM12878_1', 'PacBio_GM12878_2']
# full_df.drop(dc_datasets, inplace=True, axis=1) # drop datasets we don't want
full_df = full_df.loc[full_df[datasets].sum(axis=1) != 0] # drop transcripts with no reads in datasets we do want
full_df = full_df.loc[full_df['transcript_ID'].isin(human_tids)] # drop ebv transcripts
full_df['ebv'] = 'Human' # give human/ebv designation
full_df = full_df.loc[full_df.transcript_novelty != 'Genomic']

# drop genomic transcripts from the ebv dataset (b/c it's not pre-filtered)
ebv_df = ebv_df.loc[ebv_df.transcript_novelty != 'Genomic']
ebv_df['ebv'] = 'EBV' # human/ebv designation

# merge transcript df so TPMs can be calculated correctly
t_df = pd.concat([full_df, ebv_df])

# combine datasets
combine_datasets = True
if combine_datasets:
    t_df['combined'] = t_df['PacBio_GM12878_1']+t_df['PacBio_GM12878_2']
    datasets = ['combined']

# # make sure the concatenation worked
# print(t_df.loc[t_df['transcript_ID'] == 121].head())
# print(ebv_df.loc[ebv_df['transcript_ID'] == 121].head())
# print(full_df.loc[full_df['transcript_ID'] == 121].head())

# get tpms and number of human transcripts for 
# each dataset and for full/ebv
for d in datasets:

    # raw TPM
    tpm, t_df = get_tpm(t_df, d)

    # log2TPM
    log, t_df = get_log_tpm(t_df, d, 0)

    # sanity check - sum of all TPMs for each sample
    print('TPM total for {}: {}'.format(d, str(t_df[tpm].sum())))

    human_df = t_df.loc[(t_df[d] != 0) & (t_df['ebv'] == 'Human')]
    ebv_df = t_df.loc[(t_df[d] != 0) & (t_df['ebv'] == 'EBV')]
    n_human_transcripts = len(human_df.index)
    n_ebv_transcripts = len(ebv_df.index)
    print('Number of human transcripts in {}: {}'.format(d, str(n_human_transcripts)))
    print('Number of EBV transcripts in {}: {}'.format(d, str(n_ebv_transcripts)))

    # add columns for number of dataset human/ebv transcripts
    n_transcripts_col = 'n_'+d
    t_df[n_transcripts_col] = t_df.apply(lambda x:\
     n_human_transcripts if x['ebv'] == 'Human' else n_ebv_transcripts, axis=1)

    # add heights geom_text locations for dataset/human/ebv transcripts
    human_height = t_df.loc[t_df.ebv == 'Human'][log].max()+1
    ebv_height = t_df.loc[t_df.ebv == 'EBV'][log].max()+1
    height_col = d+'_height'
    t_df[height_col] = t_df.apply(lambda x:\
     human_height if x.ebv == 'Human' else ebv_height, axis=1)

    # print(human_height)
    # print(ebv_height)

# print(t_df.head())
# print(t_df.tail())

# write gene and transcript tables to a csv
t_df['dot_size'] = t_df.apply(lambda x: 1 if x['ebv'] == 'EBV' else 0.6, axis=1)
t_df['alpha'] = t_df.apply(lambda x: 0.5 if x['ebv'] == 'EBV' else 0.2, axis=1)


bname = get_basename(ebv_ab)
odir = format_odir(os.path.dirname(ebv_ab))
# odir = make_dated_folder(odir,bname)
to = odir+'ebv_human_transcript_abundance.csv'
transcript_outfile = open('ebv_human_transcript_abundance.csv', 'w')
t_df.to_csv(transcript_outfile, sep=',', index=False)

## get transcript tpms without filtering for bioreps to use for gene tpms
# read in the unfiltered datasets
full_df = pd.read_csv(full_unf_ab, sep='\t')
ebv_df = pd.read_csv(ebv_unf_ab, sep='\t')

# reformat human table
# dc_datasets = ['D4', 'D5', 'D10', 'D11']
datasets = ['PacBio_GM12878_1', 'PacBio_GM12878_2']
# full_df.drop(dc_datasets, inplace=True, axis=1) # drop datasets we don't want
full_df = full_df.loc[full_df['transcript_ID'].isin(human_tids)] # drop ebv transcripts
full_df['ebv'] = 'Human' # give human/ebv designation
full_df = full_df.loc[full_df.transcript_novelty != 'Genomic']

# drop genomic transcripts from the ebv dataset (b/c it's not pre-filtered)
ebv_df = ebv_df.loc[ebv_df.transcript_novelty != 'Genomic']
ebv_df['ebv'] = 'EBV' # human/ebv designation

# merge transcript df so TPMs can be calculated correctly
t_df = pd.concat([full_df, ebv_df])

# combine datasets
combine_datasets = True
if combine_datasets:
    t_df['combined'] = t_df['PacBio_GM12878_1']+t_df['PacBio_GM12878_2']
    datasets = ['combined']

# # make sure the concatenation worked
# print(t_df.loc[t_df['transcript_ID'] == 121].head())
# print(ebv_df.loc[ebv_df['transcript_ID'] == 121].head())
# print(full_df.loc[full_df['transcript_ID'] == 121].head())

# get tpms and number of human transcripts for 
# each dataset and for full/ebv
for d in datasets:

    # raw TPM
    tpm, t_df = get_tpm(t_df, d)

    # log2TPM
    log, t_df = get_log_tpm(t_df, d, 0)

    # sanity check - sum of all TPMs for each sample
    print('TPM total for {}: {}'.format(d, str(t_df[tpm].sum())))

    human_df = t_df.loc[(t_df[d] != 0) & (t_df['ebv'] == 'Human')]
    ebv_df = t_df.loc[(t_df[d] != 0) & (t_df['ebv'] == 'EBV')]
    n_human_transcripts = len(human_df.index)
    n_ebv_transcripts = len(ebv_df.index)
    print('Number of human transcripts in {}: {}'.format(d, str(n_human_transcripts)))
    print('Number of EBV transcripts in {}: {}'.format(d, str(n_ebv_transcripts)))

    # add columns for number of dataset human/ebv transcripts
    n_transcripts_col = 'n_'+d
    t_df[n_transcripts_col] = t_df.apply(lambda x:\
     n_human_transcripts if x['ebv'] == 'Human' else n_ebv_transcripts, axis=1)

    # add heights geom_text locations for dataset/human/ebv transcripts
    human_height = t_df.loc[t_df.ebv == 'Human'][log].max()+1
    ebv_height = t_df.loc[t_df.ebv == 'EBV'][log].max()+1
    height_col = d+'_height'
    t_df[height_col] = t_df.apply(lambda x:\
     human_height if x.ebv == 'Human' else ebv_height, axis=1)

# get gene tpms
cols = []
for d in datasets:
    cols.append(d)
    cols.append('TPM_'+d)
g_df = t_df.groupby(['gene_ID', 'gene_novelty', 'ebv'])[cols].sum()
g_df.reset_index(inplace=True)

# # make sure the groupby worked
# print(g_df.loc[g_df['gene_ID'] == 16].head())
# print(t_df.loc[t_df['gene_ID'] == 16].head())
# print(t_df.loc[t_df['gene_ID'] == 16].head())

# get tpms, heights, and numbers for gene
for d in datasets:
    # log2TPM
    log, g_df = get_log_tpm(g_df, d, 0)

    human_df = g_df.loc[(g_df[d] != 0) & (g_df['ebv'] == 'Human')]
    ebv_df = g_df.loc[(g_df[d] != 0) & (g_df['ebv'] == 'EBV')]
    n_human_genes = len(human_df.index)
    n_ebv_genes = len(ebv_df.index)
    print('Number of human genes in {}: {}'.format(d, str(n_human_genes)))
    print('Number of EBV genes in {}: {}'.format(d, str(n_ebv_genes)))

    # add columns for number of dataset human/ebv genes
    n_genes_col = 'n_'+d
    g_df[n_genes_col] = g_df.apply(lambda x:\
     n_human_genes if x['ebv'] == 'Human' else n_ebv_genes, axis=1)

    # add heights geom_text locations for dataset/human/ebv transcripts
    human_height = g_df.loc[g_df.ebv == 'Human'][log].max()+0.4
    ebv_height = g_df.loc[g_df.ebv == 'EBV'][log].max()+0.4
    height_col = d+'_height'
    g_df[height_col] = g_df.apply(lambda x:\
     human_height if x.ebv == 'Human' else ebv_height, axis=1)

    print(human_height)
    print(ebv_height)

print(g_df.head())
print(g_df.tail())


# add different dot sizes for human/ebv
# t_df['dot_size'] = t_df.apply(lambda x: 1 if x['ebv'] == 'EBV' else 0.6, axis=1)
g_df['dot_size'] = g_df.apply(lambda x: 1 if x['ebv'] == 'EBV' else 0.6, axis=1)

# t_df['alpha'] = t_df.apply(lambda x: 0.5 if x['ebv'] == 'EBV' else 0.2, axis=1)
g_df['alpha'] = g_df.apply(lambda x: 0.5 if x['ebv'] == 'EBV' else 0.2, axis=1)

# # rename gene/transcript novelty columns
# t_df.rename(index=str, columns={\
#     'transcript_novelty':'Isoform Type'}, inplace=True)
# g_df.rename(index=str, columns={\
#     'gene_novelty':'Gene Type'}, inplace=True)


# write gene table to a csv 

go = odir+'ebv_human_gene_abundance.csv'
gene_outfile = open('ebv_human_gene_abundance.csv', 'w')
g_df.to_csv(gene_outfile, sep=',', index=False)

# make graphs
cmd = 'Rscript plot_ebv_v_human_abundances.R --gene_csv {} --transcript_csv {}'\
	  ' --datasets {}'.format(go, to, ','.join(datasets))
# print(cmd)

