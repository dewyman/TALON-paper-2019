
#!/bin/bash
#$ -q sam128
#$ -N D4-TC-1_A01
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

set -e
module load dwyman/anaconda/3
source activate runTC
cd /pub/dwyman/TALON-paper-2019/pipeline/D4/TC_v1.0.7/1_A01
time python /data/users/dwyman/TranscriptClean-1.0.7/TranscriptClean.py --sam /pub/dwyman/TALON-paper-2019/pipeline/D4/Minimap2/1_A01/mapped_FLNC_noScaff.sam --genome /pub/dwyman/TALON-paper-2019/refs/hg38/hg38.fa --spliceJns /pub/dwyman/TALON-paper-2019/refs/TranscriptClean/gencode_v29_SJs.tsv --outprefix /pub/dwyman/TALON-paper-2019/pipeline/D4/TC_v1.0.7/1_A01/1_A01 --variants /pub/dwyman/TALON-paper-2019/refs/TranscriptClean/00-common_all.vcf.gz --primaryOnly
samtools view -H /pub/dwyman/TALON-paper-2019/pipeline/D4/TC_v1.0.7/1_A01/1_A01_clean.sam > /pub/dwyman/TALON-paper-2019/pipeline/D4/TC_v1.0.7/1_A01/canonical.sam
samtools view ${OUTDIR}/${group}/${group}_clean.sam | awk '{if($(NF-1) !~ 0) print $0}' >> /pub/dwyman/TALON-paper-2019/pipeline/D4/TC_v1.0.7/1_A01/canonical.sam
samtools view -bS /pub/dwyman/TALON-paper-2019/pipeline/D4/TC_v1.0.7/1_A01/canonical.sam | samtools sort > /pub/dwyman/TALON-paper-2019/pipeline/D4/TC_v1.0.7/1_A01/sorted_canonical.bam
samtools index /pub/dwyman/TALON-paper-2019/pipeline/D4/TC_v1.0.7/1_A01/sorted_canonical.bam
touch /pub/dwyman/TALON-paper-2019/pipeline/D4/TC_v1.0.7/1_A01/done
source deactivate
