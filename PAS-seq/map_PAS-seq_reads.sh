#!/bin/bash
#$ -q sam128
#$ -pe one-node-mpi 16
#$ -R y
#$ -N STAR_PAS-seq
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

cd /dfs2/pub/dwyman/TALON-paper-2019/PAS-seq/data
mkdir mapped_PAS

STAR --readFilesIn \
    data/PAS1_S1_L001_R1_001.fastq.gz,data/PAS1_S1_L002_R1_001.fastq.gz,data/PAS1_S1_L003_R1_001.fastq.gz,data/PAS1_S1_L004_R1_001.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix mapped/PAS_ \
    --genomeDir /dfs2/pub/dwyman/TALON-paper-2019/refs/hg38/STAR_hg38_ENCODE \
    --runThreadN 16
