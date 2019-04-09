mkdir -p data/ENCODE
cd data/ENCODE

# GM12878 CAGE IDR peaks
wget https://www.encodeproject.org/files/ENCFF853HOH/@@download/ENCFF853HOH.bed.gz
gunzip ENCFF853HOH.bed.gz
mv ENCFF853HOH.bed GM12878_CAGE.bed

# K562 CAGE IDR peaks
wget https://www.encodeproject.org/files/ENCFF698DQS/@@download/ENCFF698DQS.bed.gz
gunzip ENCFF698DQS.bed.gz
mv ENCFF698DQS.bed K562_CAGE.bed

# HepG2 CAGE IDR peaks
wget https://www.encodeproject.org/files/ENCFF246WDH/@@download/ENCFF246WDH.bed.gz
gunzip ENCFF246WDH.bed.gz
mv ENCFF246WDH.bed HepG2_CAGE.bed 
