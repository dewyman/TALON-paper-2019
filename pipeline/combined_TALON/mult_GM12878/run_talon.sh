#!/bin/bash
#$ -q sam128
#$ -N TALON
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y
set -e
module load dwyman/anaconda/3
source activate mypython3.7.2

cd /pub/dwyman/TALON-paper-2019/pipeline/combined_TALON/mult_GM12878
time python /pub/dwyman/TALON/talon.py --f config.csv \
                                       --db mult_GM12878.db \
                                       --build hg38  --cov 0.9 --identity 0 \
                                       --o mult_GM12878
source deactivate

