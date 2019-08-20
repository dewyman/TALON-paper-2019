# Extract the last exon of each GTF transcript and record in BED format

from optparse import OptionParser

def getOptions():
    parser = OptionParser()

    parser.add_option("--f", dest = "infile",
        help = "GTF file", metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def output_last_entry(o, chrom, transcript_exons, transcript_strand):
    """ Write the final exon of the provided entry to the outfile """

    if transcript_strand == None:
        return

    if transcript_strand == "+":
        last_exon = transcript_exons[-1]
    elif transcript_strand == "-":
        last_exon = transcript_exons[0]

    start = str(last_exon[0])
    end = str(last_exon[-1])

    o.write("\t".join([chrom, start, end, ".", ".", transcript_strand]) + "\n")
    return

def main():

    options = getOptions()
    infile = options.infile
    fname = options.outprefix + "_last_exons.bed"

    chrom = None
    transcript_strand = None
    transcript_exons = []

    o = open(fname, 'w')
    with open(infile, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"): continue

            entry = line.split('\t')
            entry_type = entry[2]

            if entry_type == "gene": continue

            if entry_type == "transcript":
                output_last_entry(o, chrom, transcript_exons, transcript_strand)

                # Check strand
                chrom = entry[0]
                transcript_strand = entry[6]
                transcript_exons = []

            if entry_type == "exon":
                start = int(entry[3]) - 1
                end = int(entry[4])
                transcript_exons.append([start, end])   
            
    o.close()

if __name__ == '__main__':
    main()
