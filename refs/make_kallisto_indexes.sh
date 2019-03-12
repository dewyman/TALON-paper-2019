#!/bin/bash
#$ -q sam128
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y
#$ -N kallisto_index
#$ -pe one-node-mpi 8
#$ -R y


module load kallisto/0.43.1
mkdir -p /pub/dwyman/TALON-paper-2019/refs/Kallisto
cd /pub/dwyman/TALON-paper-2019/refs/Kallisto

transcriptomeFasta=../gencode.v29.transcripts.fa.gz
indexName=gencode_v29.idx

# Run indexing step on the transcriptome
kallisto index -i $indexName $transcriptomeFasta

transcriptomeFasta=../gencode.vM20.transcripts.fa.gz
indexName=gencode_vM20.idx

# Run indexing step on the transcriptome
kallisto index -i $indexName $transcriptomeFasta
