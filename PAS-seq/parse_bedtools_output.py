# Consolidates outputs from running Bedtools Intersect on transcript ends
# against polyA sites.
# Output format: transcript_ID, PAS support (yes/no)

from optparse import OptionParser

def getOptions():
    parser = OptionParser()

    parser.add_option("--f", dest = "infile",
        help = "Bedtools intersect output", metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options


def main():
    options = getOptions()
    infile = options.infile
    outprefix = options.outprefix
    transcript_seen = {} # To avoid duplicates

    o = open(outprefix + "_PAS_results.csv", 'w')
    o.write(",".join(["transcript_ID", "PAS_support"]) + "\n")
    

    # Iterate over bedtools results
    with open(infile, 'r') as f:
        for line in f:
            line = line.strip()
            entry = line.split("\t")

            transcript_ID = entry[3]
            intersect = entry[-1]

            if transcript_ID in transcript_seen:
                continue
            if intersect != ".":
                pas_support = "yes"
            else:
                pas_support = "no"
            transcript_seen[transcript_ID] = 1

            o.write(",".join([transcript_ID, pas_support]) + "\n")

    o.close()

if __name__ == '__main__':
    main()
