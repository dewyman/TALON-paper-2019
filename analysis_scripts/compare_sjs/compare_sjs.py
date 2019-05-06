# Takes two splice junction files and determines how many junctions are shared 
# across both (or not).

from optparse import OptionParser
import os
import sys

def getOptions():
    parser = OptionParser()

    parser.add_option("--short", dest = "shortread",
        help = "Short read derived SJs", metavar = "FILE", type = "string")
    parser.add_option("--long", dest = "longread",
        help = "Long read derived SJs", metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def create_sj_set(sj_file):
    """ Given a STAR-formatted splice jn file, extract the chromosome, start,
        end, and strand of each junction. Use these to create a concatenated 
        ID, and put them into a set. """

    sj_set = set()

    with open(sj_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t') 
            
            chrom = line[0]
            start = line[1]
            end = line[2]
            strand = line[3]

            sj_entry = "_".join([chrom, start, end, strand])
            sj_set.add(sj_entry)
            
    return sj_set

def main():

    options = getOptions()
    short_read_sjs = options.shortread
    long_read_sjs = options.longread

    short_read_set = create_sj_set(short_read_sjs)
    long_read_set = create_sj_set(long_read_sjs)

    short_only = len(short_read_set - long_read_set)
    long_only = len(long_read_set - short_read_set)
    in_both = len(short_read_set & long_read_set)
    perc_long_supported = (in_both*100)/len(long_read_set)
  
    o = open(options.outprefix + "_sj_comparison.csv", 'w')
    o.write(",".join(["short_only", str(short_only)]) + "\n")
    o.write(",".join(["long_only", str(long_only)]) + "\n")
    o.write(",".join(["both", str(in_both)]) + "\n")  
    o.write(",".join(["percent_long_sjs_supported", str(perc_long_supported)]) + "\n")
    o.close()

if __name__ == '__main__':
    main()
