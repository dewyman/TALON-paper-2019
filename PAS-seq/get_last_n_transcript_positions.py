# Given a transcript bed file, create intervals of desired length from the 
# ends (ie select the last n positions).

from optparse import OptionParser
import os
import sys
import subprocess
from pathlib import Path
script_dir = os.path.dirname(os.path.realpath(__file__))
utils_dir = "/".join(script_dir.split("/")[0:-1] + ["RNA-PET"])
sys.path.append(script_dir)
sys.path.append(utils_dir)
import create_intervals as cI

def getOptions():
    parser = OptionParser()

    parser.add_option("--bed", dest = "bed",
        help = "Transcript bed file", metavar = "FILE", type = "string")
    parser.add_option("--maxdist", dest = "maxdist",
        help = "Distance from transcript end to use when creating interval",
        metavar = "FILE", type = "int")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def make_end_interval(entry, dist):
    """ Extract end position of entry and return interval of size dist
    """
    strand = entry[5]
    transcript_end = int(entry[2])

    interval_start, interval_end = cI.create_end_piece(transcript_end, strand, dist)
    return interval_start, interval_end

def main():
    options = getOptions()
    bed_file = options.bed
    dist = int(options.maxdist)
    outprefix = options.outprefix

    # Open bed files
    start_bed = open(outprefix + "_transcript_starts.bed", 'w')
    end_bed = open(outprefix + "_transcript_ends.bed", 'w')

    # Iterate over BED file
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.strip()
            entry = line.split("\t")
            chromosome = entry[0]
            name = entry[3]
            strand = entry[5]

            # Skip chrM and chrEBV
            if chromosome in ["chrM", "chrEBV"]:
                continue

            interval_start, interval_end = make_end_interval(entry, dist)

            # Write BED entry for the end
            end_entry = [ chromosome, str(interval_start),
                          str(interval_end), name, "0", strand ]
            end_bed.write("\t".join(end_entry) + "\n")

if __name__ == '__main__':
    main()
