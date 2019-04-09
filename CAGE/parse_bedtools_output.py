# Consolidates outputs from running Bedtools Intersect on transcript starts/ends
# against CAGE start and end sites.
# Output format: transcript_ID, CAGE support (yes/no)

from optparse import OptionParser

def getOptions():
    parser = OptionParser()

    parser.add_option("--f", dest = "infile",
        help = "Bedtools intersect output for the starts", metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options


def main():
    options = getOptions()
    start_file = options.infile
    outprefix = options.outprefix
    transcript_seen = {} # To avoid duplicates

    o = open(outprefix + "_CAGE_results.csv", 'w')
    o.write(",".join(["transcript_ID", "CAGE_support"]) + "\n")
    

    # Iterate over bedtools results
    with open(infile, 'r') as f:
        for line in f:
            line = line.strip()
            entry = line.split("\t")

            transcript_ID = entry[3]
            intersect = entry[-1]

            if transcript_ID in transcript_seen:
                continue
            if intersect != "-1":
                cage_support = "yes"
            else:
                cage_support = "no"
            transcript_seen[transcript_ID] = 1

            o.write(",".join([transcript_ID, cage_support]) + "\n")

    o.close()

if __name__ == '__main__':
    main()
