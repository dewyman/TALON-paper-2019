# Download common variants to use with TranscriptClean
# The subset of 00-All categorized as common (minor allele frequency >= 0.01 in at least one of 26 major populations, with at least two unrelated individuals having the minor allele)". I'm downloading this because I want to test TranscriptClean on more than on VCF file to make sure nothing breaks. https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/
mkdir -p TranscriptClean
cd TranscriptClean
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz

# Change chromosome convention so that chromosomes start with 'chr'
zcat 00-common_all.vcf.gz | \
    awk '{if($0 !~ /^#/ && $0 !~ /^chr/) print "chr"$0; else print $0}' \
    > tmp_00-common_all.vcf
gzip tmp_00-common_all.vcf
mv tmp_00-common_all.vcf.gz 00-common_all.vcf.gz
