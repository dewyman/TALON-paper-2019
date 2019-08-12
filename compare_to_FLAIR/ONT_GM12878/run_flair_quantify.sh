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
#$ -N ONT_quant

module load samtools
module load bedtools
module load freese/minimap/2-2.15
mkdir -p tmp

# cd /pub/dwyman/TALON-paper-2019/compare_to_FLAIR/GM12878

config=config.tsv
isoforms=flair.collapse.isoforms.fa

python ~/mortazavi_lab/bin/flair/flair.py quantify -i $isoforms \
                              -t 4 \
                              -r $config \
                              --tpm \
                              --temp_dir tmp
