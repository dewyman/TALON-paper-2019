# Converts a GTF file to a BED file

GTF=$1
out=$2

awk -F'\t' -v OFS='\t' '{ split($9,a,";");
          for (i=1; i <= length(a); i++)
              if(a[i] ~ "gene_id")
                  { split(a[i], b, "\""); print $1,$4-1,$5,b[2],".",$7,"."}}' $GTF \
    | bedtools sort > $out

