# Given an RNA-PET bed file, extract the start and end position of each (paying 
# attention to strand), and put the starts and ends in different bed files.

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

    parser.add_option("--rnapet", dest = "rnapet",
        help = "RNA-PET clusters in BED file format", metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def make_intervals(entry):
    """ Extract start and end position of entry and return each as bed interval 
    """
    strand = entry[5]
    rna_pet_start = int(entry[1])
    rna_pet_end = int(entry[2])

    # Set start and end based on strand
    if strand == "+":
        start_interval_1, start_interval_2 = cI.create_interval(rna_pet_start,
                                                                "left", 0)
        end_interval_1, end_interval_2 = cI.create_interval(rna_pet_end,
                                                            "right", 0)
    elif strand == "-":
        start_interval_1, start_interval_2 = cI.create_interval(rna_pet_end,
                                                                "right", 0) 
        end_interval_1, end_interval_2 = cI.create_interval(rna_pet_start,
                                                            "left", 0)
    start_tuple = (start_interval_1, start_interval_2)
    end_tuple = (end_interval_1, end_interval_2)
    return start_tuple, end_tuple

def main():
    options = getOptions()
    rna_pet_file = options.rnapet
    outprefix = options.outprefix

    # Open bed files
    start_bed = open(outprefix + "_RNA-PET_starts.bed", 'w')
    end_bed = open(outprefix + "_RNA-PET_ends.bed", 'w')

    # Iterate over BED file
    with open(rna_pet_file, 'r') as f:
        for line in f:
            line = line.strip()
            entry = line.split("\t")
            chromosome = entry[0]
            name = entry[3]
            strand = entry[5]
            
            start_interval, end_interval = make_intervals(entry)

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
