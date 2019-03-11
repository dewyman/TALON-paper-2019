set -e
mkdir -p hg38
cd hg38

# Download ENCODE human reference genome
wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz

# Remove extra information from fasta headers
awk '{print $1}' GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta > hg38.fa

# Index
module load samtools/1.3
samtools faidx hg38.fa

# Download index for STAR
wget https://www.encodeproject.org/files/ENCFF742NER/@@download/ENCFF742NER.tar.gz
tar -xvzf ENCFF742NER.tar.gz
mv out STAR_hg38_ENCODE
