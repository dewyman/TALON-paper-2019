# Given a transcript bed file, select the last n positions of the sequence and
# orient from 5' to 3'

from optparse import OptionParser
import os
import sys
import subprocess
from pathlib import Path
import pybedtools
from pyfasta import Fasta
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
        metavar = "FILE", type = "int"),
    parser.add_option("--genome", dest = "genome",
        help = "Fasta reference genome to pull sequence from",
        metavar = "FILE", type = "string"),
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


def fetch_sequence(chrom, start, end, strand, genome):
    """ Given a BED region, fetch its sequence. If it is on the minus strand,
        then reverse-complement the sequence. """

    seq = genome.sequence({'chr': chrom, 'start': start, 'stop': end, 
                      'strand': strand}, one_based=False)
    return seq 

def main():
    options = getOptions()
    bed_file = options.bed
    dist = int(options.maxdist)
    genome_file = options.genome
    outprefix = options.outprefix

    # Read genome
    genome = Fasta(genome_file)

    # Open bed files
    #end_bed = open(outprefix + "_transcript_ends.bed", 'w')

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
            seq = fetch_sequence(chromosome, interval_start, interval_end, strand, genome)
            # print(seq)
            print("AATAAA" in seq \
               or "ATTAAA" in seq \
               or "AAGAAA" in seq \
               or "AATACA" in seq \
               or "AATAGA" in seq \
               or "AATATA" in seq \
               or "AATGAA" in seq \
               or "ACTAAA" in seq \
               or "AGTAAA" in seq \
               or "CATAAA" in seq \
               or "GATAAA" in seq \
               or "TATAAA" in seq \
               or "TTTAAA" in seq)

            # Write BED entry for the end
            #end_entry = [ chromosome, str(interval_start),
            #              str(interval_end), name, "0", strand ]
            #end_bed.write("\t".join(end_entry) + "\n")

if __name__ == '__main__':
    main()
