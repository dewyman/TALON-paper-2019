#!/bin/bash
#$ -q sam
#$ -N TALON
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y
set -e
module load dwyman/anaconda/3
source activate mypython3.7.2
cd /pub/dwyman/TALON-paper-2019/pipeline/combined_TALON/
time python /pub/dwyman/TALON/talon.py --f config.csv \
                                       --db full_gencode_v29_2019-03-12.db \
                                       --build hg38  --cov 0.9 --identity 0 \
                                       --o combined
source deactivate

