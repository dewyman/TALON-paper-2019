#!/bin/bash
#$ -q sam128
#$ -pe one-node-mpi 4
#$ -R y
#$ -N flair_align_PacBio_GM12878_1
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load samtools
module load bedtools

cd /pub/dwyman/TALON-paper-2019/compare_to_FLAIR/GM12878/PacBio_GM12878_1

genome=/pub/dwyman/TALON-paper-2019/refs/hg38/hg38.fa
reads=ENCFF281TNJ.fastq

python ~/flair/flair.py align -g $genome \
                              -r $reads \
                              -t 4 \
                              -m ~/minimap2-2.15 \
                              -v1.3 

