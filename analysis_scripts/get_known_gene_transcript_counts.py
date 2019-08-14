import pandas as pd 
import sys

df = pd.read_csv(sys.argv[1], sep='\t')

cols = ['annot_gene_id', 'annot_transcript_id', 'transcript_novelty']
data_cols = df.columns.tolist()[11:]
cols.extend(df.columns.tolist()[11:])

print(data_cols)

# only keep cols we care about
df = df[cols]

# get all columns with at least one read for that gene/transcript
df = df[(df[data_cols] != 0).any(axis=1)]

# remove genomic transcripts
# df = df.loc[df.transcript_novelty != 'Genomic']

# print number of unique known genes
known_genes_df = df[df.annot_gene_id.str.contains('ENSG')]
print('known genes')
print(len(known_genes_df.annot_gene_id.unique().tolist()))

# number of unique transcripts
known_transcripts_df = df[df.annot_transcript_id.str.contains('ENST')]
print('known transcripts')
print(len(known_transcripts_df.annot_transcript_id.unique().tolist()))