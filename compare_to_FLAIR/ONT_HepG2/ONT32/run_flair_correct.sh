#!/bin/bash
#$ -q som,bio,free64
#$ -pe one-node-mpi 16
#$ -R y
#$ -m ea
#$ -cwd
#$ -j y
#$ -o /data/users/freese/mortazavi_lab/qsub_output
#$ -e /data/users/freese/mortazavi_lab/qsub_output
#$ -ckpt restart
#$ -N ONT33_correct

module load bedtools
module load samtools 
# cd /pub/dwyman/TALON-paper-2019/compare_to_FLAIR/HepG2/D8

gtf=~/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf
query=flair.aligned.bed
chromSizes=~/mortazavi_lab/ref/hg38/hg38.chrom.sizes
genome=~/mortazavi_lab/ref/hg38/hg38.fa

python ~/mortazavi_lab/bin/flair/flair.py correct -f $gtf \
                              -q $query \
                              -t 4 \
                              -c $chromSizes \
                              -g $genome
