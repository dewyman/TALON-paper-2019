# A small pipeline that compares the ends of TALON transcripts to
# starts derived from the PAS-seq assay, and plots the results by
# novelty group.

from optparse import OptionParser
import os
import sys
import subprocess
from pathlib import Path
script_dir = os.path.dirname(os.path.realpath(__file__))
utils_dir = "/".join(script_dir.split("/")[0:-1] + ["RNA-PET"])
sys.path.append(script_dir)
sys.path.append(utils_dir)

def getOptions():
    parser = OptionParser()

    parser.add_option("--gtf", dest = "gtf",
        help = "TALON GTF file", metavar = "FILE", type = "string")
    parser.add_option("--maxdist", dest = "maxdist",
        help = "Distance from transcript end in which to search for PAS motif",
        type = "int")
    parser.add_option("--genome", dest = "genome",
        help = "Reference genome for sequence",
        metavar = "FILE")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def main():
    options = getOptions()
    gtf = options.gtf
    genome = options.genome
    max_dist = options.maxdist
    outprefix = options.outprefix
    name = outprefix.split("/")[-1]

    # Step 0: set up directory structure
    try:
        outdir = "/".join(outprefix.split("/")[0:-1])
        subprocess.check_output(["mkdir", "-p", outdir + "/transcript_beds"])
    except Exception as e:
        print(e)
        sys.exit("Something went wrong while initializing dirs")

    # First step: Create a bed file from the transcript GTF along with metadata
    try:
        out = outdir + "/transcript_beds/" + name
        subprocess.check_output(["python", utils_dir +"/talon_GTF_2_transcript_bed.py",
                                 "--gtf", gtf, "--o",  out])
    except Exception as e:
        print(e)
        sys.exit("Something went wrong with talon_GTF_2_transcript_bed.py run")

    # Next, take the transcript bed file and search the end for PAS motifs.
    try:
        out = outdir + "/transcript_beds/" + name
        bedfile = outdir + "/transcript_beds/"+ name + ".bed"
        subprocess.check_output(["python", script_dir +"/check_last_n_transcript_seq_for_PA_motif.py",
                                 "--bed", bedfile, "--maxdist", str(max_dist),
                                 "--genome", genome,
                                 "--o",  outprefix])
    except Exception as e:
        print(e)
        sys.exit("Something went wrong with check_last_n_transcript_seq_for_PA_motif.py run")

if __name__ == '__main__':
    main()
