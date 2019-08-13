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
#$ -N ONT_collapse
                    
module load samtools
module load bedtools
module load freese/minimap/2-2.15

reads=ONT32-ONT33-concat.fastq
gtf=~/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf
query=ONT32-ONT33_flair_all_corrected.psl
genome=~/mortazavi_lab/ref/hg38/hg38.fa

python ~/mortazavi_lab/bin/flair/flair.py collapse -f $gtf \
                              -q $query \
                              -t 4 \
                              -g $genome \
                              -r $reads

