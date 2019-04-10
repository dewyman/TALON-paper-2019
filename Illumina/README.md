A place to put raw data files etc. pertaining to our Illumina-PacBio comparisons.

PacBio only detects transcripts that have (a) a polyA tail, and (b) are over 300 bp (unless you change the pipeline defaults). Because of this, we filter the Illumina gene set to consist of genes that it is possible for PacBio to detect. We use the following procedure:  
* Start with Kallisto file  
* Remove transcripts with a TPM < 1
* Remove transcripts with length < 300 bp  
* Remove mitochondrial genes (lack a polyA tail):  
```
mitochondrial_blacklist <- c("MT-TF", "MT-RNR1", "MT-TV", "MT-RNR2", "MT-TL1", "MT-ND1", "MT-TI", "MT-TQ", "MT-TM", "MT-ND2", "MT-TW", "MT-TA", "MT-TN", "MT-TC", "MT-TY", "MT-CO1", "MT-TS1", "MT-TD", "MT-CO2", "MT-TK", "MT-ATP8", "MT-ATP6", "MT-CO3", "MT-TG", "MT-ND3", "MT-TR", "MT-ND4L", "MT-ND4", "MT-TH", "MT-TS2", "MT-TL2", "MT-ND5", "MT-ND6", "MT-CYB", "MTATP6P1")
```

After performing these filtering steps, we should normalize the expression values so that they sum to 1M again.

When performing an analysis focusing on genes rather than transcripts, we change the filtering slightly: we aggregate the transcripts for each gene before applying the TPM cutoff. The reason: we know that short read quantification estimates on the transcript level are approximate, so we want to consider all of the reads mapping to the gene when doing this filtering.

Scripts for these things:
```
TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_transcripts.R
TALON-paper-2019/analysis_scripts/filter_kallisto_illumina_genes.R
```
