#!/bin/bash
#$ -q sam128
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y
#$ -N kallisto_GM12878
#$ -pe one-node-mpi 16
#$ -R y


module load kallisto/0.43.1
cd /dfs2/pub/dwyman/TALON-paper-2019/Illumina/GM12878
mkdir -p Kallisto

illumina1=GM12878_rep1_illumina_1.fastq.gz
illumina2=GM12878_rep1_illumina_2.fastq.gz
transcriptomeFasta=/dfs2/pub/dwyman/TALON-paper-2019/refs/gencode.v29.transcripts.fa.gz
indexName=/dfs2/pub/dwyman/TALON-paper-2019/refs/Kallisto/gencode_v29.idx

kallisto quant -i $indexName -t 16 -o Kallisto -b 100 $illumina1 $illumina2
