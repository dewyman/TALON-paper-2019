from optparse import OptionParser
import os

def getOptions():
    parser = OptionParser()

    parser.add_option("--f", dest = "infile",
        help = "Input GTF file of genomic transcripts", metavar = "FILE", type = "string")
    parser.add_option("--p", dest = "prefix", help = "Prefix for intermediate files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def main():

    options = getOptions()
    infile = options.infile
    bed = options.prefix + "_genomic.bed"

    # Get genomic transcripts in BED format
    cmd = """awk -v OFS='\t' '{if($3 == "transcript") print $1,$4-1,$5,".",".",$7}' """ + infile + " > " + bed
    os.system(cmd)

    # Bedtools intersect it
    btools_out = options.prefix + "_nGenomic_intersect_lastExons.bed"
    bedtools_cmd = """bedtools intersect -a GM12878_pacbio_repro_genomic_talon.bed \
                      -b gencode_v29_last_exons.bed \
                      -u \
                      -s | wc -l > """ + btools_out
    os.system(bedtools_cmd)

    # Now, collect the results for output
    with open(btools_out) as f:
        overlap_size = int(f.readline().strip())

    total_genomic = sum(1 for line in open(bed))
    percent_overlap = round(overlap_size*100./total_genomic)

    print("\t".join([str(overlap_size), str(total_genomic), str(percent_overlap) + "%"]))

if __name__ == '__main__':
    main()
