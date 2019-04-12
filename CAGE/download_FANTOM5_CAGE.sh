mkdir -p data/FANTOM5
cd data/ENCODE

# hg19 robust CAGE peaks
wget http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_coord.bed.gz
gunzip hg19.cage_peak_phase1and2combined_coord.bed.gz
mv hg19.cage_peak_phase1and2combined_coord.bed hg19_CAGE.bed
