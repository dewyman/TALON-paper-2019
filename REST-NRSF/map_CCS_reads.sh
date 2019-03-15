
#!/bin/bash
#$ -q sam128
#$ -pe one-node-mpi 16
#$ -R y
#$ -N Minimap2
#$ -M dwyman@uci.edu
#$ -m ea
#$ -cwd
#$ -j y

module load samtools/1.3

cd /pub/dwyman/TALON-paper-2019/REST-NRSF 

reads_A=/share/samdata/dwyman/PB69/CCS/1_A01/ccs.bam
reads_B=/share/samdata/dwyman/PB69/CCS/2_B01/ccs.bam

# Convert reads to fasta and concatenate
fasta=HepG2_D4_ccs_reads.fasta
samtools fasta -n $reads_A  > $fasta
samtools fasta -n $reads_B  >> $fasta

# Run Minimap2 on the reads
mapped_ccs=HepG2_D4_ccs_reads.sam
/data/users/dwyman/minimap2-2.15/minimap2 -t 16 -ax splice -uf --secondary=no -C5 \
           /pub/dwyman/TALON-paper-2019/refs/hg38/hg38.fa       \
           $fasta \
           > $mapped_ccs         \
           2> mapped_ccs.log
samtools view -bS $mapped_ccs > /pub/public-www/dwyman/HepG2_D4_ccs_reads.bam
samtools index /pub/public-www/dwyman/HepG2_D4_ccs_reads.bam
