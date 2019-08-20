# Genomic transcripts
The goal here is to characterize the properties of the 'genomic' transcript category.

```
wget http://crick.bio.uci.edu/talon_supplement/tables/190615_supplemental_tables_submission/S27_full_gencode_v29_pb_ont_talon_abundance.tsv
```

1. Create a list of genomic transcripts that were reproducible in the GM12878 PacBio transcriptome
```
awk -v OFS=',' '{if($10 == "Genomic" && $14*$15 > 0) print $1,$2,$10}' S27_full_gencode_v29_pb_ont_talon_abundance.tsv > GM12878_PacBio_reproducible_genomic.csv
```

2. Use this file and the corresponding TALON database to generate a GTF file of only reproducible 
genomic transcripts.
```
python /pub/dwyman/TALON/post-TALON_tools/create_GTF_from_database.py \
          --db full_gencode_v29_2019-05-24.db \
          -a gencode_v29 \
          -b hg38 \
          --whitelist GM12878_PacBio_reproducible_genomic.csv \
          --datasets GM12878_pacbio_datasets.txt \
          --o GM12878_pacbio_repro_genomic
```

3. Generate antisense mapping file needed for the CAGE/PAS/RNA-PET plotting

## Map antisense TALON genes to their sense counterparts
```
python /pub/dwyman/TALON/post-TALON_tools/map_antisense_genes_to_sense.py \
    --db full_gencode_v29_2019-05-24.db \
    --a gencode_v29 \
    --o tier1
```


4. Perform FANTOM CAGE analysis
```
source activate mypython3.7.2
mkdir -p GM12878-PacBio/CAGE
python ../CAGE/run_CAGE_analysis.py \
        --gtf GM12878_pacbio_repro_genomic_talon.gtf \
        --cage ../CAGE/data/FANTOM5/hg38_CAGE.bed \
        --maxdist 100 \
        --o GM12878-PacBio/CAGE/GM12878-genomic

Rscript ../analysis_scripts/plot_support_by_novelty_type.R \
    --f GM12878-PacBio/CAGE/GM12878-genomic_CAGE_results.csv \
    --t CAGE \
    --novelty GM12878-PacBio/CAGE/transcript_beds/GM12878-genomic_novelty.csv \
    --abundance S27_full_gencode_v29_pb_ont_talon_abundance.tsv \
    --d1 PacBio_GM12878_1 --d2 PacBio_GM12878_2 --as tier1_antisense_mapping.csv \
    -o GM12878-PacBio/CAGE/GM12878-genomic
```

5. Perform computational PAS analysis
```
source activate mypython3.7.2
mkdir -p GM12878-PacBio/PAS-comp
python ../PAS-computational/run_computational_PAS_analysis.py \
        --gtf GM12878_pacbio_repro_genomic_talon.gtf \
        --genome ../refs/hg38/hg38.fa \
        --maxdist 35 \
        --o GM12878-PacBio/PAS-comp/GM12878-genomic

Rscript ../analysis_scripts/plot_support_by_novelty_type.R \
    --f GM12878-PacBio/PAS-comp/GM12878-genomic_polyA_motif.csv \
    --t PAS-comp \
    --novelty GM12878-PacBio/PAS-comp/transcript_beds/GM12878-genomic_novelty.csv \
    --abundance S27_full_gencode_v29_pb_ont_talon_abundance.tsv \
    --d1 PacBio_GM12878_1 --d2 PacBio_GM12878_2 --as tier1_antisense_mapping.csv \
    -o GM12878-PacBio/PAS-comp/GM12878-genomic
```

6. Perform RNA-PET analysis
```
mkdir -p GM12878-PacBio/RNA-PET
module load bedtools
python ../RNA-PET/run_RNA-PET_analysis.py \
    --gtf GM12878_pacbio_repro_genomic_talon.gtf \
    --rnapet ../RNA-PET/data/GM12878_hg38.bed \
    --maxdist 100 \
    --o GM12878-PacBio/RNA-PET/GM12878-genomic

Rscript ../analysis_scripts/plot_support_by_novelty_type.R  \
     --f GM12878-PacBio/RNA-PET/GM12878-genomic_RNA-PET_results.csv   \
     --t RNA-PET \
     --novelty GM12878-PacBio/RNA-PET/transcript_beds/GM12878-genomic_novelty.csv \
     --abundance S27_full_gencode_v29_pb_ont_talon_abundance.tsv \
     --d1 PacBio_GM12878_1 --d2 PacBio_GM12878_2 \
     --as tier1_antisense_mapping.csv \
     -o GM12878-PacBio/RNA-PET/GM12878-genomic

```

7. Find out how many genomic transcripts overlap with the last exon of a transcript
```
# Get the last exon of every known GENCODE transcript
python extract_last_exons.py --f ../refs/gencode.v29.annotation.gtf --o gencode_v29

# Get transcripts in BED format
awk -v OFS='\t' '{if($3 == "transcript") print $1,$4-1,$5,".",".",$7}' GM12878_pacbio_repro_genomic_talon.gtf > GM12878_pacbio_repro_genomic_talon.bed

# Bedtools intersect it
source activate mypython3.7.2
bedtools intersect -a GM12878_pacbio_repro_genomic_talon.bed \
                   -b gencode_v29_last_exons.bed \
                   -u \
                   -s | wc -l 
```
