## Analysis of unknown genes detected in long reads from Tier 1 ENCODE cell lines

We are interested in what proportion of the genes belong to the following categories:
1. Multiexonic gene that overlaps a known gene, but is antisense
2. Overlaps a known gene on the same strand, but does not share any splice junctions with it
3. Intergenic
4. Monoexonic gene that overlaps a known gene, but is antisense
5. Overlaps an ALU or LINE1 repeat
6. < 1000 bp upstream of known gene (either strand)
7. < 1000 bp downstream of known gene (either strand)
8. On the EBV chromosome (which is not included in the GTF annotation)

In cases where the query overlaps a known gene, that will be prioritized over containing a repeat.

## The Data
First, we generate a file of novel transcript/gene IDs for the novel models that were pairwise reproducible in at least one of the Tier 1 cell lines. Then we use this list to create a transcriptome GTF from the PacBio data. This is the file entitled 'novel_models_talon.gtf'.
```
awk -v OFS=',' '{if($9 != "Known" && ($12*$13 > 0 || $14*$15 >0 || $16*$17 >0)) print $1,$2,$10}' ${tier1_abd_filt} > novel_models_talon.csv

talon_create_GTF \
          --db ${human_db} \
          -a gencode_v29 \
          -b hg38 \
          --whitelist novel_models_talon.csv \
          --datasets tier1_datasets.txt \
          --o novel_models
```
We also generate a file indicating whether each of these genes is monoexonic or multiexonic using the abundance table:
```
# The first sorting command sorts by column 1 first and then by column 2 numerically and in reverse order
# The second will remove duplicates taking into account just the first column
awk -v OFS='\t' '{if($9 != "Known" && ($12*$13 > 0 || $14*$15 >0 || $16*$17 >0)) print $3,$7}' ${tier1_abd_filt} \
    | sort -k1,1 -k2,2rn \
    | sort -uk1,1 \        
    | awk -v OFS='\t' '{if($2 == 1) print $1,"monoexonic"; else print $1,"multiexonic"}' \
    > n_exons_for_novel_genes.tsv
```
We also need the full transcriptome GTF, ../refs/gencode.v29.annotation.gtf.

Alu and LINE1 repeat files for the human genome (Dec. 2013 GRCh38/hg38) were downloaded in the BED format on 10-20-18 from the UCSC genome table browser as follows:
1) Navigate to the UCSC table browser and select the following options:
Repeats, RepeatMasker, rmask, BED output format.
2) After downloading the file, run the following command to extract only Alu repeats:
```
grep "Alu" hg38_repeats_rmsk.bed > hg38_alu_repeats_rmsk.bed
```
3) Run the following command to extract only LINE1 repeats:
```
zcat hg38_repeats_rmsk.bed.gz | awk '{if($4 ~ "^L1") print $0}' > hg38_LINE1_repeats_rmsk.bed
```

## Data Processing
Data processing was performed using the script run_novel_gene_analysis.sh.
To summarize, the first step is to extract the novel genes from the GTF and format them as a BED file. We must also extract the known genes in BED format so that we can compare the novel genes to them. What follows is a series of Bedtools Intersect and Bedtools Closest commands in order to compare the novel genes to the known ones. The genes are also intersected with the alu and LINE1 repeats. Then, bash commands are used to parse the Bedtools output and combine the results of different comparisons into a single file (entitled 'tier1/master_file.tsv'). 
```
source activate mypython3.7.2
# In case of Bedtools segfaults, log directly into sam128 node
./run_novel_gene_analysis.sh novel_models_talon.gtf \
                             ../refs/gencode.v29.annotation.gtf \
                             n_exons_for_novel_genes.tsv \
                             hg38_alu_repeats_rmsk.bed \
                             hg38_LINE1_repeats_rmsk.bed \
                             .
```

## Data Visualization
Plot a pie chart that illustrates what percentage of the novel genes belong to each category. It takes the file 'master_file.tsv' as input. 
```
module load R/3.5.1
Rscript plot_gene_categorization.R --f master_file.tsv --o plots
```


