# Download ENCODE human reference genome
mkdir -p mm10
cd mm10
wget https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz
gunzip mm10_no_alt_analysis_set_ENCODE.fasta.gz

# Index
module load samtools/1.3
samtools faidx mm10_no_alt_analysis_set_ENCODE.fasta

# Download index for STAR
wget https://www.encodeproject.org/files/ENCFF483PAE/@@download/ENCFF483PAE.tar.gz
tar -xvzf ENCFF483PAE.tar.gz
mv out STAR_mm10_ENCODE
