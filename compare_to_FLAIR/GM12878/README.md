# Running FLAIR on GM12878 PacBio data

FLAIR was cloned from https://github.com/BrooksLabUCSC/flair on 8/5/2019.

1. Run align and correct steps separately on replicates
```bash
qsub PacBio_GM12878_1/./run_FLAIR_align.sh
qsub PacBio_GM12878_2/./run_FLAIR_align.sh
```
```bash
qsub PacBio_GM12878_1/./run_FLAIR_correct.sh
qsub PacBio_GM12878_2/./run_FLAIR_correct.sh
```
2. Then, run collapse step on concatenated files from both reps.
```bash
cat PacBio_GM12878_1/flair_all_corrected.psl PacBio_GM12878_2/flair_all_corrected.psl > PacBio_GM12878_1-PacBio_GM12878_2_flair_all_corrected.psl
cat PacBio_GM12878_1/ENCFF281TNJ.fastq PacBio_GM12878_2/ENCFF475ORL.fastq > PacBio_GM12878_1-PacBio_GM12878_2-concat.fastq
qsub ./run_flair_collapse.sh
```
3. Finally, run quantify step. To do this, you need to create a tab-delimited config file with fields dataset name, condition, batch, and fastq reads file. This is what the GM12878 file looks like:
```bash
GM12878_Rep1	GM12878	batch1	/pub/dwyman/TALON-paper-2019/compare_to_FLAIR/PacBio_GM12878_1/ENCFF281TNJ.fastq
GM12878_Rep2	GM12878	batch1	/pub/dwyman/TALON-paper-2019/compare_to_FLAIR/PacBio_GM12878_2/ENCFF475ORL.fastq
```
```bash
qsub ./run_flair_quantify.sh
```

4. In order to determine how well PacBio + FLAIR detects genes known to be expressed in short-read data, we converted the FLAIR output to a TALON-like format, and then ran a custom R script:
```bash
python ../format_flair_matrix_like_talon.py counts_matrix.tsv counts_matrix_talon_abd.tsv

Rscript ../plot_detection_by_TPM_for_datasets.R \
      --f counts_matrix_talon_abd.tsv \
      --datasets GM12878_Rep1_GM12878_batch1,GM12878_Rep2_GM12878_batch1 \
      --ik1 ../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
      --ik2 ../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
      --color blue \
      -o .

Rscript ../pacbio_v_illumina_edgeR.R \
    --f counts_matrix_talon_abd.tsv \
    --datasets GM12878_Rep1_GM12878_batch1,GM12878_Rep2_GM12878_batch1 \
    --ik1 ../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
    --ik2 ../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
    --color blue \
    -o FLAIR
```

5. We also want to see how reproducible our datasets are as characterized by FLAIR. Run the following gene and transcript correlations to see:
```bash
Rscript ../plot_pacbio_gene_expression_corr.R \
    --f counts_matrix_talon_abd.tsv \
    --color blue \
    --d1 GM12878_Rep1_GM12878_batch1 \
    --d2 GM12878_Rep2_GM12878_batch1 \
    --celltype GM12878 \
    --d1_label "PacBio Rep1" \
    --d2_label "PacBio Rep2" \
    -o FLAIR 

Rscript ../plot_pacbio_transcript_expression_corr.R \
   --f counts_matrix_talon_abd.tsv \
   --color blue \
    --d1 GM12878_Rep1_GM12878_batch1 \
    --d2 GM12878_Rep2_GM12878_batch1 \
    --celltype GM12878 \
    --d1_label "PacBio Rep1" \
    --d2_label "PacBio Rep2" \
    --outdir FLAIR 
```

6. Finally, let's see how the gene detection sensitivity compares between FLAIR and TALON on the same datasets. We will need the PacBio GM12878 unfiltered TALON abundance file 
```bash
sup_tables=/data/users/freese/TALON_data/revisions_10-19/human_TALON/analysis/supplementary_tables/
abundance=${sup_tables}S3_GM12878_talon_abundance.tsv

Rscript ../compare_TALON_FLAIR_detection_to_Illumina.R \
    --talon ${abundance} \
    --flair counts_matrix_talon_abd.tsv \
    --talonD PacBio_GM12878_1,PacBio_GM12878_2 \
    --flairD GM12878_Rep1_GM12878_batch1,GM12878_Rep2_GM12878_batch1 \
    --ik1 ../../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
    --ik2 ../../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
    -o . 
```
