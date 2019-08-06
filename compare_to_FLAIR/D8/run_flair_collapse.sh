#!/bin/bash
#$ -q sam
#$ -pe one-node-mpi 4
#$ -R y
#$ -N D8_flair_collapse
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load samtools
module load bedtools

cd /pub/dwyman/TALON-paper-2019/compare_to_FLAIR/D8

reads=ENCFF281TNJ.fastq
gtf=/pub/dwyman/TALON-paper-2019/refs/gencode.v29.annotation.gtf
query=flair_all_corrected.psl
genome=/pub/dwyman/TALON-paper-2019/refs/hg38/hg38.fa

python ~/flair/flair.py collapse -f $gtf \
                              -q $query \
                              -t 4 \
                              -g $genome \
                              -m ~/minimap2-2.15 \
                              -r $reads

