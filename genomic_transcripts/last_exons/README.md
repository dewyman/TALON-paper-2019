# What proportion of reproducible genomic transcripts overlap with the last exon of a known transcript?

```
database=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_10-19/human_TALON/full_gencode_v29_2019-06-19.db
sup_tables=/share/crsp/lab/seyedam/share/TALON_paper_data/revisions_10-19/human_TALON/analysis/supplementary_tables
tier1_filtered=${sup_tables}/S28_full_gencode_v29_pb_ont_talon_abundance.tsv
```

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

Column numbers:  
12. PacBio_HepG2_1
13. PacBio_HepG2_2
14. PacBio_GM12878_1
15. PacBio_GM12878_2
16. PacBio_GM12878_3
17. PacBio_GM12878_4
18. PacBio_K562_1
19. PacBio_K562_2
20. ONT_HepG2_1
21. ONT_HepG2_2
22. ONT_HepG2_3
23. ONT_GM12878_1
24. ONT_GM12878_2
25. ONT_GM12878_3
26. ONT_K562_1
27. ONT_K562_2

```
awk -v OFS=',' '{if($10 == "Genomic" && $14*$15 > 0) print $1,$2,$10}' ${tier1_filtered} > PacBio_GM12878/reproducible_genomic.csv

awk -v OFS=',' '{if($10 == "Genomic" && $18*$19 > 0) print $1,$2,$10}' ${tier1_filtered}  > PacBio_K562/reproducible_genomic.csv

awk -v OFS=',' '{if($10 == "Genomic" && $12*$13 > 0) print $1,$2,$10}' ${tier1_filtered} > PacBio_HepG2/reproducible_genomic.csv

awk -v OFS=',' '{if($10 == "Genomic" && $23*$24 > 0) print $1,$2,$10}' ${tier1_filtered} > ONT_GM12878/reproducible_genomic.csv

awk -v OFS=',' '{if($10 == "Genomic" && $26*$27 > 0) print $1,$2,$10}' ${tier1_filtered} > ONT_K562/reproducible_genomic.csv

awk -v OFS=',' '{if($10 == "Genomic" && $20*$21 > 0) print $1,$2,$10}' ${tier1_filtered} > ONT_HepG2/reproducible_genomic.csv
```

4. Make dataset files
```
echo "PacBio_GM12878_1" > PacBio_GM12878/datasets.txt
echo "PacBio_GM12878_2" >> PacBio_GM12878/datasets.txt

echo "PacBio_K562_1" > PacBio_K562/datasets.txt
echo "PacBio_K562_2" >> PacBio_K562/datasets.txt

echo "PacBio_HepG2_1" > PacBio_HepG2/datasets.txt
echo "PacBio_HepG2_2" >> PacBio_HepG2/datasets.txt

echo "ONT_GM12878_1" > ONT_GM12878/datasets.txt
echo "ONT_GM12878_2" >> ONT_GM12878/datasets.txt

echo "ONT_K562_1" > ONT_K562/datasets.txt
echo "ONT_K562_2" >> ONT_K562/datasets.txt

echo "ONT_HepG2_1" > ONT_HepG2/datasets.txt
echo "ONT_HepG2_2" >> ONT_HepG2/datasets.txt
```


5. Use TALON database to make GTF for each set of genomic transcipts
```
talon_create_GTF \
          --db ${database} \
          -a gencode_v29 \
          -b hg38 \
          --whitelist PacBio_GM12878/reproducible_genomic.csv \
          --datasets PacBio_GM12878/datasets.txt \
          --o PacBio_GM12878/reproducible_genomic.gtf

talon_create_GTF \
          --db ${database} \
          -a gencode_v29 \
          -b hg38 \
          --whitelist PacBio_K562/reproducible_genomic.csv \
          --datasets PacBio_K562/datasets.txt \
          --o PacBio_K562/reproducible_genomic.gtf

talon_create_GTF \
          --db ${database} \
          -a gencode_v29 \
          -b hg38 \
          --whitelist PacBio_HepG2/reproducible_genomic.csv \
          --datasets PacBio_HepG2/datasets.txt \
          --o PacBio_HepG2/reproducible_genomic.gtf

talon_create_GTF \
          --db ${database} \
          -a gencode_v29 \
          -b hg38 \
          --whitelist ONT_GM12878/reproducible_genomic.csv \
          --datasets ONT_GM12878/datasets.txt \
          --o ONT_GM12878/reproducible_genomic.gtf

talon_create_GTF \
          --db ${database} \
          -a gencode_v29 \
          -b hg38 \
          --whitelist ONT_K562/reproducible_genomic.csv \
          --datasets ONT_K562/datasets.txt \
          --o ONT_K562/reproducible_genomic.gtf

talon_create_GTF \
          --db ${database} \
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
