# PAS-seq
We compared the 3' ends of our long-read transcripts to data from the PAS-seq assay in order to validate them.

## Data
K562 PAS-seq data was provided by the Spitale lab.
```
./download_PAS-seq.sh
```

## Analysis
Reads were mapped using STAR
```
./map_PAS-seq_reads.sh

module load samtools
samtools view data/mapped_PAS/PAS_Aligned.out.sam | awk '{print $5}' > data/mapped_PAS/MAPQ_scores.csv
Rscript plot_MAPQ_distribution.R --f data/mapped_PAS/MAPQ_scores.csv --o data/mapped_PAS/MAPQ_scores.png
```
The range of MAPQ scores in the alignments were examined. According to the STAR manual,
* 255 indicates a uniquely mapped read  
* int(-10*log10(1- 1/Nmap)) is used for multimappers  
Plotting the ditribution of MAPQ scores in the reads, we obtained the following:
![MAPQ_scores.png](./data/mapped_PAS/MAPQ_scores.png)
Based on this distribution, we elected to apply a MAPQ threshold of 255 to filter the reads (just keep the unique mappers).
```
samtools view -Sq 255 data/mapped_PAS/PAS_Aligned.out.sam > data/mapped_PAS/filt_PAS_Aligned.out.sam
```
