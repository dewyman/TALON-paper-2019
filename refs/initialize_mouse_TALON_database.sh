#!/bin/bash
#$ -q sam128
#$ -N init_mouse_db_full
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

DATE=$(date +%F\(%a\) | awk -F"(" '{print $1}')

module load dwyman/anaconda/3
source activate mypython3.7.2
cd /pub/dwyman/TALON-paper-2019/refs
mkdir -p TALON
cd TALON

time python /pub/dwyman/TALON/initialize_talon_database.py \
    --f ../gencode.vM21.annotation.gtf \
    --a gencode_vM21 \
    --g mm10 \
    --l 300 \
    --idprefix ENCODEM \
    --5p 500 \
    --3p 300 \
    --o unmodified_full_gencode_vM21_${DATE}

source deactivate
