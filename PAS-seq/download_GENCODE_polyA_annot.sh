wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.metadata.PolyA_feature.gz
gunzip gencode.v29.metadata.PolyA_feature.gz

# Convert to BED file. I'm assuming that the original file is 1-based because GENCODE's other fiels are.
# Sample line:
# ENST00000484859.1	4859	4860	chr1	141474	141475	-	polyA_site
awk -v OFS="\t" '{if($8 == "polyA_signal") print $4,$5-1,$6,$1"_"$8,".",$7}' gencode.v29.metadata.PolyA_feature > gencode.v29.metadata.PolyA_feature.bed
