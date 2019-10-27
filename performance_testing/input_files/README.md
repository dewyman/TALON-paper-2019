# Input files

To create input files for the runtime performance tests, follow these instructions:

0. Download and install seqtk as described [here](https://github.com/lh3/seqtk) 

1. Download FLNC human PacBio reads from the ENCODE portal
```
wget https://www.encodeproject.org/files/ENCFF281TNJ/@@download/ENCFF281TNJ.fastq.gz
wget https://www.encodeproject.org/files/ENCFF475ORL/@@download/ENCFF475ORL.fastq.gz
wget https://www.encodeproject.org/files/ENCFF329AYV/@@download/ENCFF329AYV.fastq.gz
wget https://www.encodeproject.org/files/ENCFF902UIT/@@download/ENCFF902UIT.fastq.gz
wget https://www.encodeproject.org/files/ENCFF763VZC/@@download/ENCFF763VZC.fastq.gz
wget https://www.encodeproject.org/files/ENCFF694INI/@@download/ENCFF694INI.fastq.gz
wget https://www.encodeproject.org/files/ENCFF427JDY/@@download/ENCFF427JDY.fastq.gz
wget https://www.encodeproject.org/files/ENCFF589SMB/@@download/ENCFF589SMB.fastq.gz
cat *.fastq.gz > combined.fq.gz
```

2. Concatenate the reads

3. Randomly select 1000 reads from the combined file. Then repeat process for 10,000, 100,000, 1M, 2.5M, and 5M reads. The -s parameter refers to the random seed.
```
SEQTK_PATH=~/seqtk

$SEQTK_PATH/seqtk sample -s100 combined.fq.gz 1000 > 1000.fq
$SEQTK_PATH/seqtk sample -s100 combined.fq.gz 10000 > 10000.fq
$SEQTK_PATH/seqtk sample -s100 combined.fq.gz 100000 > 100000.fq
$SEQTK_PATH/seqtk sample -s100 combined.fq.gz 1000000 > 1000000.fq
$SEQTK_PATH/seqtk sample -s100 combined.fq.gz 2500000 > 2500000.fq
$SEQTK_PATH/seqtk sample -s100 combined.fq.gz 5000000 > 5000000.fq
```
