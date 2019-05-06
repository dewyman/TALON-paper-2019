#!/bin/bash
#$ -q sam128
#$ -pe one-node-mpi 16
#$ -R y
#$ -N STAR-GM12878
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

sample1_reads1=/pub/dwyman/TALON-paper-2019/Illumina/GM12878/GM12878_rep1_illumina_1.fastq.gz
sample1_reads2=/pub/dwyman/TALON-paper-2019/Illumina/GM12878/GM12878_rep1_illumina_2.fastq.gz
sample2_reads1=/pub/dwyman/TALON-paper-2019/Illumina/GM12878/GM12878_rep2_illumina_1.fastq.gz
sample2_reads2=/pub/dwyman/TALON-paper-2019/Illumina/GM12878/GM12878_rep2_illumina_2.fastq.gz

module load STAR/2.5.2a
STAR --runThreadN 16 --genomeDir /pub/dwyman/TALON-paper-2019/refs/hg38/STAR_hg38_ENCODE \
     --readFilesIn $sample1_reads1,$sample1_reads2,$sample2_reads1,$sample2_reads2 \
     --sjdbGTFfile /pub/dwyman/TALON-paper-2019/refs/gencode.v29.annotation.gtf \
     --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
     --outFileNamePrefix /pub/dwyman/TALON-paper-2019/Illumina/GM12878/STAR \
     --outSAMattributes NH HI NM MD jM jI --outSAMtype SAM --readFilesCommand zcat
