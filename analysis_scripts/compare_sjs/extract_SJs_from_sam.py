# A script to take all primary-mapped reads in a SAM file and report the splice
# junctions present, including the splice motif. The output format is the same
# as STAR uses:
# 1) chromosome
# 2) first intron base (1-based)
# 3) last intron base (1-based)
# 4) strand 
#       0 = undefined
#       1 = + 
#       2 = -
# 5) intron motif:
#       0 = noncanonical (not one of the following)
#       1 = GT/AG
#       2 = CT/AC
#       3 = GC/AG
#       4 = CT/GC
#       5 = AT/AC
#       6 = GT/AT
# 6) Annotation status (0/1)
# 7) Number of uniquely mapped reads on this junction
# 8) Number of multimapping reads on this junction (we set to NA)
# 9) Max spliced alignment overhang (we set to NA)

import pyfasta as pf
from optparse import OptionParser
import re

def getOptions():
    parser = OptionParser()

    parser.add_option("--sam", dest = "sam",
        help = "SAM file to extract junctions from", metavar = "FILE", type = "string")
    parser.add_option("--genome", dest = "genome",
        help = "Reference genome fasta file", metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output file",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def split_cigar(cigar):
    """ Takes CIGAR string from SAM and splits it into two lists:
        one with capital letters (match operators), and one with
        the number of bases that each operation applies to. """

    alignTypes = re.sub('[0-9]', " ", cigar).split()
    counts = re.sub('[A-Z]', " ", cigar).split()
    counts = [int(i) for i in counts]

    return alignTypes, counts

def get_introns(start, cigar):
    """ Compute position of introns. This is done by stepping through the CIGAR
        string, where introns are represented by the N operation. 
        Note that positions refer to start and endpoints of introns, not exons,
        so adjustments are needed to avoid an off-by-one error if you want exons.
        Second note: coordinates will be in ascending order, not necessarily in 
        their sense order.
    """
    operations, counts = split_cigar(cigar)
    intron_coords = []
    genomePos = start

    # Iterate over cigar operations
    for op,ct in zip(operations, counts):
        if op == "N":
            # This is an intron
            intronStart = genomePos
            intronEnd = genomePos + ct - 1

            intron_coords.append(intronStart)
            intron_coords.append(intronEnd)

        if op not in ["S", "I"]:
            genomePos += ct

    return intron_coords

def get_SJs_from_read(read_fields, genome):
    """ Working from the SAM fields, first compute intron positions and splice 
        motifs (1-based). Then, iterate over pairs of positions to create the 
        splice junction coords. 
        Returns: A list of tuples
                    [ (chr, start, end, strand, splice_motif), .... ]
                 One tuple per junction.
    """

    sjs = []
    
    # Compute each splice junction position and splice motif
    FLAG = int(read_fields[1])
    chrom = read_fields[2]
    start = int(read_fields[3])
    if FLAG == 0:
        strand = 1
    elif FLAG == 16:
        strand = 2
    else:
        strand = 0

    CIGAR = read_fields[5]
    intron_coords = get_introns(start, CIGAR)
    print(intron_coords)

    
    

def main():
    options = getOptions()
    sam_file = options.sam
    fasta = options.genome
    outprefix = options.outprefix

    # Read reference genome
    genome = pf.Fasta(fasta)
    
    # Initialize data structures

    # Iterate over each transcript. Disregard if it isn't a primary mapper
    with open(sam_file, 'r') as f:
        for read in f:
            read = read.strip() 
            if read.startswith("@"): continue

            read_fields = read.split('\t')
            FLAG = int(read_fields[1])
            CHROM = read_fields[2]
            if (FLAG != 0 and FLAG != 16) or CHROM == "*": continue

            # Compute each splice junction position and splice motif
            list_of_sjs = get_SJs_from_read(read_fields, genome)

            # Add to dict to track number of reads containing the junction


    # Output all of the junctions to a file     



if __name__ == '__main__':
    main() 
