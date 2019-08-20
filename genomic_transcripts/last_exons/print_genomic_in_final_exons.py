from optparse import OptionParser
import os

def getOptions():
    parser = OptionParser()

    parser.add_option("--f", dest = "infile",
        help = "Input GTF file of genomic transcripts", metavar = "FILE", type = "string")
    parser.add_option("--e", dest = "exons",
        help = "Input BED file of final exons", metavar = "FILE", type = "string")
    parser.add_option("--p", dest = "prefix", help = "Prefix for intermediate files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def main():

    options = getOptions()
    infile = options.infile
    exons = options.exons
    bed = options.prefix + "genomic.bed"

    # Get genomic transcripts in BED format
    cmd = """awk -v OFS='\t' '{if($3 == "transcript") print $1,$4-1,$5,".",".",$7}' """ + infile + " > " + bed
    os.system(cmd)

    # Bedtools intersect it
    btools_out = options.prefix + "nGenomic_intersect_lastExons.bed"
    bedtools_cmd = """bedtools intersect -a %s \
                      -b %s \
                      -u \
                      -s | wc -l > %s""" % (bed, exons, btools_out)
    os.system(bedtools_cmd)

    # Now, collect the results for output
    with open(btools_out) as f:
        overlap_size = int(f.readline().strip())

    total_genomic = sum(1 for line in open(bed))
    percent_overlap = round(overlap_size*100./total_genomic)

    print("\t".join([options.prefix, str(overlap_size), str(total_genomic), str(percent_overlap) + "%"]))

if __name__ == '__main__':
    main()
