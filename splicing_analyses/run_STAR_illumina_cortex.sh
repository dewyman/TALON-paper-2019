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
#$ -N cortex_st

DPATH=~/mortazavi_lab/bin/TALON-paper-2019/analysis_scripts/compare_sjs/mouse_brain/
curpath=~/mortazavi_lab/bin/TALON-paper-2019/splicing_analyses/technology_SJ_comparison/
reads="${DPATH}Cx1A_R1.fastq.gz,${DPATH}Cx2A_R1.fastq.gz ${DPATH}Cx1A_R2.fastq.gz,${DPATH}Cx2A_R2.fastq.gz"

module load STAR/2.5.2a
STAR \
	--runThreadN 8 \
	--genomeDir ~/mortazavi_lab/ref/hg38/STAR_hg38_ENCODE/ \
	--readFilesIn $reads \
	--sjdbGTFfile ~/mortazavi_lab/ref/gencode.v29/gencode.v29.annotation.gtf \
	--outFilterType BySJout \
	--outFilterMultimapNmax 20 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverLmax 0.04 \
	--alignIntronMin 20 \
	--alignIntronMax 1000000 \
	--alignMatesGapMax 1000000 \
	--outFileNamePrefix ${curpath}cortex_aligned \
	--outSAMattributes NH HI NM MD jM jI \
	--outSAMtype SAM \
	--readFilesCommand zcat
