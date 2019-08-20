# What proportion of reproducible genomic transcripts overlap with the last exon of a known transcript?

1. Get the last exon of every known GENCODE transcript
```
python extract_last_exons.py --f ../../refs/gencode.v29.annotation.gtf --o gencode_v29
```

2. Set up dirs
```
mkdir -p PacBio_GM12878
mkdir -p PacBio_K562
mkdir -p PacBio_HepG2
mkdir -p ONT_GM12878
mkdir -p ONT_K562
mkdir -p ONT_HepG2
```

3. Create a list of genomic transcripts that were reproducible in the:
    - GM12878 PacBio transcriptome (R1 & R2)
    - K562 PacBio transcriptome (R1 & R2)
    - HepG2 PacBio transcriptome (R1 & R2)
    - GM12878 ONT transcriptome (R1 & R2)
    - K562 ONT transcriptome (R1 & R2)
    - HepG2 ONT transcriptome (R1 & R3)

```
awk -v OFS=',' '{if($10 == "Genomic" && $14*$15 > 0) print $1,$2,$10}' ../S27_full_gencode_v29_pb_ont_talon_abundance.tsv > PacBio_GM12878/reproducible_genomic.csv

awk -v OFS=',' '{if($10 == "Genomic" && $16*$17 > 0) print $1,$2,$10}' ../S27_full_gencode_v29_pb_ont_talon_abundance.tsv > PacBio_K562/reproducible_genomic.csv

awk -v OFS=',' '{if($10 == "Genomic" && $12*$13 > 0) print $1,$2,$10}' ../S27_full_gencode_v29_pb_ont_talon_abundance.tsv > PacBio_HepG2/reproducible_genomic.csv

awk -v OFS=',' '{if($10 == "Genomic" && $19*$20 > 0) print $1,$2,$10}' ../S27_full_gencode_v29_pb_ont_talon_abundance.tsv > ONT_GM12878/reproducible_genomic.csv

awk -v OFS=',' '{if($10 == "Genomic" && $24*$25 > 0) print $1,$2,$10}' ../S27_full_gencode_v29_pb_ont_talon_abundance.tsv > ONT_K562/reproducible_genomic.csv

awk -v OFS=',' '{if($10 == "Genomic" && $22*$23 > 0) print $1,$2,$10}' ../S27_full_gencode_v29_pb_ont_talon_abundance.tsv > ONT_HepG2/reproducible_genomic.csv
```

4. Make dataset files
```
echo "D8" > PacBio_GM12878/datasets.txt
echo "D9" >> PacBio_GM12878/datasets.txt

echo "D10" > PacBio_K562/datasets.txt
echo "D11" >> PacBio_K562/datasets.txt

echo "D4" > PacBio_HepG2/datasets.txt
echo "D5" >> PacBio_HepG2/datasets.txt

echo "ONT25" > ONT_GM12878/datasets.txt
echo "ONT24" >> ONT_GM12878/datasets.txt

echo "ONT18" > ONT_K562/datasets.txt
echo "ONT31" >> ONT_K562/datasets.txt

echo "ONT32" > ONT_HepG2/datasets.txt
echo "ONT33" >> ONT_HepG2/datasets.txt
```


5. Use TALON database to make GTF for each set of genomic transcipts
```
python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist PacBio_GM12878/reproducible_genomic.csv \
          --datasets PacBio_GM12878/datasets.txt \
          --o PacBio_GM12878/reproducible_genomic.gtf

python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist PacBio_K562/reproducible_genomic.csv \
          --datasets PacBio_K562/datasets.txt \
          --o PacBio_K562/reproducible_genomic.gtf

python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist PacBio_HepG2/reproducible_genomic.csv \
          --datasets PacBio_HepG2/datasets.txt \
          --o PacBio_HepG2/reproducible_genomic.gtf

python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist ONT_GM12878/reproducible_genomic.csv \
          --datasets ONT_GM12878/datasets.txt \
          --o ONT_GM12878/reproducible_genomic.gtf

python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist ONT_K562/reproducible_genomic.csv \
          --datasets ONT_K562/datasets.txt \
          --o ONT_K562/reproducible_genomic.gtf

python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db ../full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist ONT_HepG2/reproducible_genomic.csv \
          --datasets ONT_HepG2/datasets.txt \
          --o ONT_HepG2/reproducible_genomic.gtf
```


6. Run script to intersect each set of genomic transcripts with the final exon set from GENCODE and record the results.
```
python print_genomic_in_final_exons.py --f PacBio_GM12878/reproducible_genomic_talon.gtf --e gencode_v29_last_exons.bed --p PacBio_GM12878/ > genomic_transcripts_in_final_exons.tsv

python print_genomic_in_final_exons.py --f PacBio_K562/reproducible_genomic_talon.gtf --e gencode_v29_last_exons.bed --p PacBio_K562/ >> genomic_transcripts_in_final_exons.tsv

python print_genomic_in_final_exons.py --f PacBio_HepG2/reproducible_genomic_talon.gtf --e gencode_v29_last_exons.bed --p PacBio_HepG2/ >> genomic_transcripts_in_final_exons.tsv

python print_genomic_in_final_exons.py --f ONT_GM12878/reproducible_genomic_talon.gtf --e gencode_v29_last_exons.bed --p ONT_GM12878/ >> genomic_transcripts_in_final_exons.tsv

python print_genomic_in_final_exons.py --f ONT_K562/reproducible_genomic_talon.gtf --e gencode_v29_last_exons.bed --p ONT_K562/ >> genomic_transcripts_in_final_exons.tsv

python print_genomic_in_final_exons.py --f ONT_HepG2/reproducible_genomic_talon.gtf --e gencode_v29_last_exons.bed --p ONT_HepG2/ >> genomic_transcripts_in_final_exons.tsv
```
