#!/bin/bash
#$ -q sam
#$ -pe one-node-mpi 4
#$ -R y
#$ -N D9_flair_correct
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load bedtools

cd /pub/dwyman/TALON-paper-2019/compare_to_FLAIR/GM12878/D9

gtf=/pub/dwyman/TALON-paper-2019/refs/gencode.v29.annotation.gtf
query=flair.aligned.bed
chromSizes=/pub/dwyman/TALON-paper-2019/refs/hg38/hg38.chrom.sizes
genome=/pub/dwyman/TALON-paper-2019/refs/hg38/hg38.fa

python ~/flair/flair.py correct -f $gtf \
                              -q $query \
                              -t 4 \
                              -c $chromSizes \
                              -g $genome

