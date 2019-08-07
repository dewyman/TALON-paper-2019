


In order to determine how well PacBio + FLAIR detects genes known to be expressed in short-read data, we converted the FLAIR output to a TALON-like format, and then ran a custom R script:
```
python format_flair_matrix_like_talon.py counts_matrix.tsv counts_matrix_talon_abd.tsv

Rscript plot_detection_by_TPM_for_datasets.R \
      --f counts_matrix_talon_abd.tsv \
      --datasets GM12878_Rep1_GM12878_batch1,GM12878_Rep2_GM12878_batch1 \
      --ik1 ../Illumina/GM12878/Kallisto/Rep1/abundance.tsv \
      --ik2 ../Illumina/GM12878/Kallisto/Rep2/abundance.tsv \
      --color blue \
      -o .
```
