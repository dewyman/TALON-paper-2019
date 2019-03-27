# Download RNA-PET data from the ENCODE portal. Mapped to hg19

mkdir -p data
cd data

# HepG2: clone
wget https://www.encodeproject.org/files/ENCFF001TIR/@@download/ENCFF001TIR.bed.gz
mv ENCFF001TIR.bed.gz HepG2.bed.gz

# GM12878: clone-free
wget https://www.encodeproject.org/files/ENCFF001TIL/@@download/ENCFF001TIL.bed.gz
mv ENCFF001TIL.bed.gz GM12878.bed.gz

# K562: clone-free
wget https://www.encodeproject.org/files/ENCFF001TJA/@@download/ENCFF001TJA.bed.gz
mv ENCFF001TJA.bed.gz K562.bed.gz
