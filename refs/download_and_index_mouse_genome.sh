# Download ENCODE human reference genome
# This is time consuming- run in screen

mkdir -p mm10
cd mm10
wget https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz
gunzip mm10_no_alt_analysis_set_ENCODE.fasta.gz

# Remove extra information from fasta headers
awk '{print $1}' mm10_no_alt_analysis_set_ENCODE.fasta > mm10.fa

# Index
module load samtools/1.3
samtools faidx mm10.fa

# Download index for STAR
wget https://www.encodeproject.org/files/ENCFF483PAE/@@download/ENCFF483PAE.tar.gz
tar -xvzf ENCFF483PAE.tar.gz
mv out STAR_mm10_ENCODE
