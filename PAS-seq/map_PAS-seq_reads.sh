#!/bin/bash
#$ -q sam128
#$ -pe one-node-mpi 32
#$ -R y
#$ -N STAR_PAS-seq
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load STAR/2.6.0c
cd /dfs3/pub/dwyman/TALON-paper-2019/PAS-seq/data
mkdir -p mapped_PAS

STAR --readFilesIn \
    data/PAS1_S1_L001_R1_001.fastq.gz,data/PAS1_S1_L002_R1_001.fastq.gz,data/PAS1_S1_L003_R1_001.fastq.gz,data/PAS1_S1_L004_R1_001.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix mapped_PAS/PAS_ \
    --genomeDir /dfs3/pub/dwyman/TALON-paper-2019/refs/hg38/STAR_hg38_ENCODE \
    --runThreadN 32
