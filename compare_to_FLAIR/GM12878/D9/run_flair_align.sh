#!/bin/bash
#$ -q sam
#$ -pe one-node-mpi 4
#$ -R y
#$ -N flair_align_D9
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load samtools
module load bedtools

cd /pub/dwyman/TALON-paper-2019/compare_to_FLAIR/GM12878/D9

genome=/pub/dwyman/TALON-paper-2019/refs/hg38/hg38.fa
reads=ENCFF475ORL.fastq

python ~/flair/flair.py align -g $genome \
                              -r $reads \
                              -t 4 \
                              -m ~/minimap2-2.15 \
                              -v1.3 

