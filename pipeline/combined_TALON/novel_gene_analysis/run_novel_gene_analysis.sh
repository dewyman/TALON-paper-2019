
module load bedtools/2.25.0

# Command line options
novel_GTF=$1 # GTF file of novel genes
known_GTF=$2 # GTF file of known genes
ALU_FILE=$3 # BED file of Alu repeats
LINE1_FILE=$4 # BED file of LINE1 repeats
PROJDIR=$5 # Directory in which to conduct the analysis (i.e. write output files)

# Setup
WORKDIR=$PROJDIR/intermediate_files
plot_dir=$PROJDIR/plots
mkdir -p $plot_dir
mkdir -p $WORKDIR


# ------- Extract known and novel genes in BED format from the GTF -----------
known_genes=$WORKDIR/known_genes.bed
novel_genes=$WORKDIR/novel_genes.bed

# Select lines containing known genes, and format into BED. Extract the TALON 
# gene ID to serve as the name field
awk '{if ($3 == "gene") print $0}' $known_GTF \
    | grep 'gene_status "KNOWN"' | ./gtf_2_bed.sh - $known_genes

# Do the same for novel genes
awk '{if ($3 == "gene") print $0}' $novel_GTF \
    | grep 'gene_status "NOVEL"' | ./gtf_2_bed.sh - $novel_genes

# Determine whether each novel gene is multiexonic or monoexonic
n_exon_file=$WORKDIR/gene_multiexonic.tsv
paste -d '\t' <(zcat $novel_GTF | awk '{if ($3 == "exon") print $0}' \
                  | grep 'gene_status "NOVEL"' \
                  | awk -F'\t' -v OFS='\t' '{ split($9,a,";");
                  for (i=1; i <= length(a); i++)
                      if(a[i] ~ "talon_gene")
                          { split(a[i], b, "\""); print b[2]}}') \
              <(zcat $novel_GTF | awk '{if ($3 == "exon") print $0}' \
                  | grep 'gene_status "NOVEL"' \
                  | awk -F'\t' -v OFS='\t' '{ split($9,a,";");
                  for (i=1; i <= length(a); i++)
                      if(a[i] ~ "exon_number")
                          { split(a[i], b, "\""); print b[2]}}') \
           | sort -k1,1 -k2,2rn | sort -uk1,1 \
           | awk '{if ($2 > 1) print $1,"multiexonic"; else print $1,"monoexonic"}' \
           > $n_exon_file

# -------------------------- Intersections -----------------------------------
overlap_ss=$WORKDIR/overlap_sameStrand_gene.txt
overlap_as=$WORKDIR/overlap_antisense_gene.txt
closest_ss_gene=$WORKDIR/closest_sameStrand_gene_dist.txt
closest_as_gene=$WORKDIR/closest_antisense_gene_dist.txt
alu_overlap=$WORKDIR/alu_overlap.txt
line1_overlap=$WORKDIR/line1_overlap.txt

# For each novel gene, find the greatest amount of overlap with a known gene
# on the same strand. When using the intersect command, it is necessary to remove
# duplicates in cases where there was more than one match, keeping the biggest overlap
bedtools intersect -a $novel_genes -b $known_genes -s -wao \
    | awk '{print $4,$NF}' | sort -k1,1 -k2,2rn | sort -uk1,1 > $overlap_ss

# Greatest amount of antisense overlap
bedtools intersect -a $novel_genes -b $known_genes -S -wao \
    | awk '{print $4,$NF}' | sort -k1,1 -k2,2rn | sort -uk1,1 > $overlap_as

# For each novel gene, find the closest known gene on the same strand
# (d = 0 means overlapping)
bedtools closest -a $novel_genes -b $known_genes -s -D a -t first \
    | awk '{print $4,$NF}' > $closest_ss_gene

# Find closest antisense gene (d = 0 means overlapping)
bedtools closest -a $novel_genes -b $known_genes -S -D a -t first \
    | awk '{print $4,$NF}' > $closest_as_gene

# Find amount of overlap with alu repeats
bedtools intersect -a $novel_genes -b $ALU_FILE -s -wao \
    | awk '{print $4,$NF}' | sort -k1,1 -k2,2rn | sort -uk1,1 > $alu_overlap

# Find amount of overlap with LINE1 repeats
bedtools intersect -a $novel_genes -b $LINE1_FILE -s -wao \
    | awk '{print $4,$NF}' | sort -k1,1 -k2,2rn | sort -uk1,1 > $line1_overlap

#-------------- Join together the intersection results on gene ID --------------
master_file=$PROJDIR/master_file.tsv
column_names=("gene_ID" "chromosome" "start" "end" "known_gene_overlap" "known_gene_as_overlap" "closest_gene_dist" "closest_gene_as_dist" "alu_overlap" "line1_overlap" "n_exons")
( IFS=$'\t'; echo "${column_names[*]}" ) > $master_file

join <(awk '{print $4,$1,$2,$3}' $novel_genes | sort ) $overlap_ss \
    | join - <(sort $overlap_as) \
    | join - <(sort $closest_ss_gene) \
    | join - <(sort $closest_as_gene) \
    | join - <(sort $alu_overlap) \
    | join - <(sort $line1_overlap) \
    | join - <(sort $n_exon_file) \
    | awk -v OFS='\t' '$1=$1' >> $master_file


