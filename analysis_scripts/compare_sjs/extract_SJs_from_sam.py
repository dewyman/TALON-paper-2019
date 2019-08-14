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

def get_splice_seq(chrom, pos, bound, genome):
        """ The splice motif consists of the first two or last two bases of the
            intron (first two if bound == 0 and last two if bound == 1).
            Bound = 0 refers to the left side of the jn, and bound = 1 refers to
            the right side (forward orientation) """

        if bound == 0:
            motif = genome.sequence({'chr': chrom, 'start': pos, 
                                     'stop': pos + 1}, one_based=True)
        else:
            motif = genome.sequence({'chr': chrom, 'start': pos - 1, 
                                     'stop': pos}, one_based=True)
        return motif

def getSJMotifCode(startBases, endBases):
    """ Determines which STAR-style splice junction code applies to a splice motif """       

    motif = (startBases + endBases).upper()

    if motif == "GTAG":
        return 1
    elif motif == "CTAC":
        return 2
    elif motif == "GCAG":
        return 3
    elif motif == "CTGC":
        return 4
    elif motif == "ATAC":
        return 5
    elif motif == "GTAT":
        return 6
    else:
        return 0

def fetch_splice_motif_code(chrom, start_pos, end_pos, strand, genome):
    """ Use Pyfasta to extract the splice motif sequence based on the start
        and end of the splice junction. Then, convert this sequence motif 
        to a numeric code. """

    start_motif = get_splice_seq(chrom, start_pos, 0, genome)
    end_motif = get_splice_seq(chrom, end_pos, 1, genome)

    motif_code = getSJMotifCode(start_motif, end_motif)
    return motif_code

def create_sj_tuples(chrom, strand, intron_coords, genome):
    """ Walk through intron coord pairs to assemble SJ tuples for each:
            (chr, start, end, strand, splice_motif)
        Since we want the smallest coordinate first in each case, we do not 
        need to orient based on strand."""

    sj_tuples = []

    start_index = 0
    while start_index < len(intron_coords) -1 :
        end_index = start_index + 1

        start_pos = intron_coords[start_index]
        end_pos = intron_coords[end_index]

        motif = fetch_splice_motif_code(chrom, start_pos, end_pos, strand, genome)

        # Assemble the tuple
        sj = (chrom, start_pos, end_pos, strand, motif)
        sj_tuples.append(sj)

        start_index += 2

    return sj_tuples

def get_SJs_from_read(read_fields, genome):
    """ Working from the SAM fields, first compute intron positions and splice 
        motifs (1-based). Then, iterate over pairs of positions to create the 
        splice junction coords. 
        Returns: A list of tuples
                    [ (chr, start, end, strand, splice_motif), .... ]
                 One tuple per junction.
    """

    # Compute each splice junction position and splice motif
    FLAG = int(read_fields[1])
    chrom = read_fields[2]
    start = int(read_fields[3])
    CIGAR = read_fields[5]
    intron_coords = get_introns(start, CIGAR)

    if FLAG == 0:
        strand = 1
    elif FLAG == 16:
        strand = 2
    else:
        strand = 0 
    
    sjs = create_sj_tuples(chrom, strand, intron_coords, genome)
    return sjs

def main():
    options = getOptions()
    sam_file = options.sam
    fasta = options.genome
    outprefix = options.outprefix

    # Read reference genome
    genome = pf.Fasta(fasta)
    
    # Initialize data structures
    sj_count = {}

    # Iterate over each transcript. Disregard if it isn't a primary mapper
    print("Parsing junctions...")
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

            # Add each junction to dict in order to track number of reads 
            # that have it
            for sj in list_of_sjs:
                try:
                    sj_count[sj] += 1
                except:
                    sj_count[sj] = 1 

    # Output all of the junctions to a file
    print("Writing junctions to output...")
    o = open(outprefix + "_SJs.txt", 'w')     
    for junction, count in sj_count.items():
        cols = [str(x) for x in junction] + ["0", str(count), "NA", "NA"]
        o.write("\t".join(cols) + "\n")
    o.close()    


if __name__ == '__main__':
    main() 
