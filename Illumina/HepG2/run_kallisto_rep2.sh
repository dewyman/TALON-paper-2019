#!/bin/bash
#$ -q sam128
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y
#$ -N kallisto_HepG2_R2
#$ -pe one-node-mpi 8
#$ -R y


module load kallisto/0.43.1
cd /pub/dwyman/TALON-paper-2019/Illumina/HepG2
mkdir -p Kallisto/Rep2

illumina1=HepG2_rep2_illumina_1.fastq.gz
illumina2=HepG2_rep2_illumina_2.fastq.gz
transcriptomeFasta=/pub/dwyman/TALON-paper-2019/refs/gencode.v29.transcripts.fa.gz
indexName=/pub/dwyman/TALON-paper-2019/refs/Kallisto/gencode_v29.idx

kallisto quant -i $indexName -t 8 -o Kallisto/Rep2 -b 100 $illumina1 $illumina2
