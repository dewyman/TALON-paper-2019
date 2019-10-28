# Given a transcript bed file, create intervals of desired size around the 
# starts and ends that go in separate files. Note: the interval size refers
# to the distance on one side (ie 100 bp would mean 100 bp upstream OR 100 bp
# downstream of the initial position).

from optparse import OptionParser
import os
import sys
import subprocess
from pathlib import Path
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, os.pardir)))
import create_intervals as cI

def getOptions():
    parser = OptionParser()

    parser.add_option("--bed", dest = "bed",
        help = "Transcript bed file", metavar = "FILE", type = "string")
    parser.add_option("--maxdist", dest = "maxdist",
        help = "Distance from transcript site to use when creating interval",
        metavar = "FILE", type = "int")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def make_intervals(entry, dist):
    """ Extract start and end position of entry and return each as bed interval
    """
    strand = entry[5]
    transcript_start = int(entry[1])
    transcript_end = int(entry[2])

    # Set start and end based on strand
    if strand == "+":
        start_interval_1, start_interval_2 = cI.create_interval(transcript_start,
                                                                "left", dist)
        end_interval_1, end_interval_2 = cI.create_interval(transcript_end,
                                                            "right", dist)
    elif strand == "-":
        start_interval_1, start_interval_2 = cI.create_interval(transcript_end,
                                                                "right", dist)
        end_interval_1, end_interval_2 = cI.create_interval(transcript_start,
                                                            "left", dist)
    start_tuple = (start_interval_1, start_interval_2)
    end_tuple = (end_interval_1, end_interval_2)
    return start_tuple, end_tuple

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

            # Skip chrM
            if chromosome in ["chrM"]:
                continue

            start_interval, end_interval = make_intervals(entry, dist)

            # Write BED entry for the start
            start_entry = [ chromosome, str(start_interval[0]),
                            str(start_interval[1]), name, "0", strand ]
            start_bed.write("\t".join(start_entry) + "\n")

            # Write BED entry for the end
            end_entry = [ chromosome, str(end_interval[0]),
                          str(end_interval[1]), name, "0", strand ]
            end_bed.write("\t".join(end_entry) + "\n")

if __name__ == '__main__':
    main()
