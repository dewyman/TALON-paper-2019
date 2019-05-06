#!/bin/bash
#$ -q sam128
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y
#$ -N kallisto_K562_R1
#$ -pe one-node-mpi 16
#$ -R y


module load kallisto/0.43.1
cd /pub/dwyman/TALON-paper-2019/Illumina/K562
mkdir -p Kallisto/Rep1

illumina1=K562_rep1_illumina_1.fastq.gz
illumina2=K562_rep1_illumina_2.fastq.gz
transcriptomeFasta=/pub/dwyman/TALON-paper-2019/refs/gencode.v29.transcripts.fa.gz
indexName=/pub/dwyman/TALON-paper-2019/refs/Kallisto/gencode_v29.idx

kallisto quant -i $indexName -t 16 -o Kallisto/Rep1 -b 100 $illumina1 $illumina2
