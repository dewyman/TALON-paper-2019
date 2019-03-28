# Download RNA-PET data from the ENCODE portal. Mapped to hg19
# Modify file to remove the score field because it crashes liftOver

mkdir -p data
cd data

# HepG2: clone
wget https://www.encodeproject.org/files/ENCFF001TIR/@@download/ENCFF001TIR.bed.gz
gunzip ENCFF001TIR.bed.gz 
awk -v OFS="\t" '{print $1,$2,$3,$4,0,$6}' ENCFF001TIR.bed > HepG2.bed
rm ENCFF001TIR.bed

# GM12878: clone-free
wget https://www.encodeproject.org/files/ENCFF001TIL/@@download/ENCFF001TIL.bed.gz
gunzip ENCFF001TIL.bed.gz
awk -v OFS="\t" '{print $1,$2,$3,$4,0,$6}' ENCFF001TIL.bed > GM12878.bed
rm ENCFF001TIL.bed

# K562: clone-free
wget https://www.encodeproject.org/files/ENCFF001TJA/@@download/ENCFF001TJA.bed.gz
gunzip ENCFF001TJA.bed.gz
awk -v OFS="\t" '{print $1,$2,$3,$4,0,$6}' ENCFF001TJA.bed > K562.bed
rm ENCFF001TJA.bed
