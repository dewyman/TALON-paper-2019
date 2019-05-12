# Directions for running an isoform switching analysis with IsoformSwitchAnalyzeR

## Installs

1. Install Bioconductor in an R session. I am using v3.5.1 of R
```
if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    BiocManager::install()
}
```
2. Now install IsoformSwitchAnalyzeR in an R session
```
BiocManager::install("IsoformSwitchAnalyzeR")
```

## File prep
You will need the abundance file from the TALON run, along with the GTF.
You also need a fasta file for your GTF annotation. You can get that using a TransCoder utility
```
module load freese/TransDecoder
${UTILPATH}gtf_genome_to_cdna_fasta.pl $GTF $REF > ${OPATH}${BNAME}.fasta
```
Here, $REF refers to the reference genome fasta file.


