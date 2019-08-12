i# Running FLAIR on K562 PacBio data

FLAIR was cloned from https://github.com/BrooksLabUCSC/flair on 8/5/2019.

1. Run align and correct steps separately on replicates
```
qsub ONT18/./run_FLAIR_align.sh
qsub ONT31/./run_FLAIR_align.sh
```
```
qsub ONT18/./run_FLAIR_correct.sh
qsub ONT31/./run_FLAIR_correct.sh
```
2. Then, run collapse step on concatenated files from both reps.
```
cat ONT18/flair_all_corrected.psl ONT31/flair_all_corrected.psl > ONT18-ONT31_flair_all_corrected.psl
cat ONT18/ONT_K562_1.fastq ONT31/ONT_K562_2.fastq > ONT18-ONT31-concat.fastq
qsub ./run_flair_collapse.sh
```
3. Finally, run quantify step. To do this, you need to create a tab-delimited config file with fields dataset name, condition, batch, and fastq reads file. This is what the K562 file looks like:
```
K562_Rep1	K562	batch1	/pub/dwyman/TALON-paper-2019/compare_to_FLAIR/ONT18/ONT_K562_3.fastq
K562_Rep2	K562	batch1	/pub/dwyman/TALON-paper-2019/compare_to_FLAIR/ONT31/ONT_K562_1.fastq
```
```
qsub ./run_flair_quantify.sh
```

4. In order to determine how well PacBio + FLAIR detects genes known to be expressed in short-read data, we converted the FLAIR output to a TALON-like format, and then ran a custom R script:
```
python ../format_flair_matrix_like_talon.py counts_matrix.tsv counts_matrix_talon_abd.tsv

Rscript ../plot_detection_by_TPM_for_datasets.R \
      --f counts_matrix_talon_abd.tsv \
      --datasets K562_ONT_Rep1_K562_batch1,K562_ONT_Rep2_K562_batch1 \
      --ik1 ~/mortazavi_lab/bin/TALON-paper-2019/Illumina/K562/Kallisto/Rep1/abundance.tsv \
      --ik2 ~/mortazavi_lab/bin/TALON-paper-2019/Illumina/K562/Kallisto/Rep2/abundance.tsv \
      --color red \
      -o FLAIR
`
Rscript ../pacbio_v_illumina_edgeR.R \
    --f counts_matrix_talon_abd.tsv \
    --datasets K562_ONT_Rep1_K562_batch1,K562_ONT_Rep2_K562_batch1 \
    --ik1 ~/mortazavi_lab/bin/TALON-paper-2019/Illumina/K562/Kallisto/Rep1/abundance.tsv \
    --ik2 ~/mortazavi_lab/bin/TALON-paper-2019/Illumina/K562/Kallisto/Rep2/abundance.tsv \
    --color red \
    -o FLAIR
```
<img align="left" width="500" src="FLAIR/gene_detection_by_TPM.png">
<img align="left" width="500" src="FLAIR/edgeR_pacbio_illumina_gene_MA_plot.png">

