# A small pipeline that compares the starts and ends of TALON transcripts to
# starts and ends derived from the RNA-PET assay, and plots the results by 
# novelty group.

from optparse import OptionParser
import os
import sys
import subprocess
from pathlib import Path
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, os.pardir)))

def getOptions():
    parser = OptionParser()

    parser.add_option("--gtf", dest = "gtf",
        help = "TALON GTF file", metavar = "FILE", type = "string")
    parser.add_option("--rnapet", dest = "rnapet",
        help = "RNA-PET clusters in BED file format", metavar = "FILE", type = "string")
    parser.add_option("--maxdist", dest = "maxdist",
        help = "Distance (+ or -) from transcript site in which to search for RNA-PET hit", 
        metavar = "FILE", type = "int")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def main():
    options = getOptions()
    gtf = options.gtf
    rna_pet_file = options.rnapet
    max_dist = options.maxdist
    outprefix = options.outprefix

    # First step: Create a bed file from the transcript GTF along with metadata
    try:
        subprocess.check_output(["python", "talon_GTF_2_transcript_bed.py",
                                 "--gtf", gtf, "--o",  outprefix])
    except Exception as e:    
        print(e)
        sys.exit("Something went wrong with talon_GTF_2_transcript_bed.py run")

    # Next, take the transcript bed file and create intervals of specified size
    # around the starts and ends. These go into separate files.
    try:
        bedfile = outprefix + ".bed"
        subprocess.check_output(["python", "get_transcript_start_end_intervals.py",
                                 "--bed", bedfile, "--maxdist", str(max_dist),
                                 "--o",  outprefix])
    except Exception as e:
        print(e)
        sys.exit("Something went wrong with get_transcript_start_end_intervals.py run")

    # Extract the start and end points (len 1) of each RNA-PET cluster. Put them
    # into separate files
    try:
        subprocess.check_output(["python", "get_RNA_PET_starts_and_ends.py",
                                 "--rnapet", rna_pet_file, "--o",  outprefix])
    except Exception as e:
        print(e)
        sys.exit("Something went wrong with get_RNA_PET_starts_and_ends.py run")

    exit()
    # Run Bedtools intersect on (a) transcript starts and RNA-PET starts
    #                           (b) transcript ends and RNA-PET ends
    try:
        pass
    except Exception as e:
        pass

    # Take the Bedtools output and determine whether each transcript start-end 
    # pair matched with at least one RNA-PET start-end pair or not. 
    try:
        pass
    except Exception as e:
        pass

    # Generate an RNA-PET coverage plot based on the output file from the 
    # previous step




if __name__ == '__main__':
    main()
