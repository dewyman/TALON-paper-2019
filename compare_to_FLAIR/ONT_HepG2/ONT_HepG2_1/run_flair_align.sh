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
#$ -N ONT_HepG2_1_align


module load samtools
module load bedtools
module load freese/minimap/2-2.15

# cd /pub/dwyman/TALON-paper-2019/compare_to_FLAIR/HepG2/PacBio_GM12878_1

genome=~/mortazavi_lab/ref/hg38/hg38.fa
reads=ONT_HepG2_1.fastq

python ~/mortazavi_lab/bin/flair/flair.py align -g $genome \
                              -r $reads \
                              -t 4 \
                              -v1.3
