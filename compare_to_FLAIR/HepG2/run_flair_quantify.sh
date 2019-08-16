#!/bin/bash
#$ -q sam
#$ -pe one-node-mpi 4
#$ -R y
#$ -N HepG2_flair_quantify
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load samtools
module load bedtools
mkdir -p tmp

cd /pub/dwyman/TALON-paper-2019/compare_to_FLAIR/HepG2

config=config.tsv
isoforms=flair.collapse.isoforms.fa

python ~/flair/flair.py quantify -i $isoforms \
                              -t 4 \
                              -m ~/minimap2-2.15 \
                              -r $config \
                              --tpm \
                              --temp_dir tmp
