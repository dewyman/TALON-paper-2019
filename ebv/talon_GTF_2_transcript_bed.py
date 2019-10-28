# Take a GTF file output from TALON and create a BED file containing 
# the transcripts from the GTF. Also create a metadata file that
# tracks the types of novelty that the transcripts have.

from optparse import OptionParser

def getOptions():
    parser = OptionParser()

    parser.add_option("--gtf", dest = "gtf",
        help = "GTF file", metavar = "FILE", type = "string")

    parser.add_option("--o", dest = "outprefix", help = "Prefix for output files",
        metavar = "FILE", type = "string")


    (options, args) = parser.parse_args()
    return options

def create_BED_entry(gtf_transcript):
    """ Given a GTF transcript (in list form), create a BED entry. This entails:
        1. Convert coordinates from 1-based to 0-based
        2. Extract unique transcript identifier to use in name field
        3. Extract other attributes (ie chromosome, strand) """

    chromosome = gtf_transcript[0]
    start_1b = int(gtf_transcript[3])
    end_1b = int(gtf_transcript[4])
    strand = gtf_transcript[6]
    meta = gtf_transcript[-1] 
   
    # Convert coords to 0-based
    start_0b = start_1b - 1
    end_0b = end_1b

    # Extract unique transcript identifier
    # Parse description
    transcript_ID = parse_out_transcript_ID(meta)

    # assemble BED
    bed = [ chromosome, start_0b, end_0b, transcript_ID, "0", strand ]
    return bed

def parse_out_transcript_ID(gtf_metadata):
    """ Parse GTF metadata in order to extract the transcript ID """
    transcript_ID = (gtf_metadata.split('transcript_id "')[1]).split('";')[0]
    return(transcript_ID)

def create_metadata_entry(gtf_transcript):
    """ Given a GTF transcript (in list form), determine which novelty types
        the transcript has. """

    meta = gtf_transcript[-1]
    transcript_ID = parse_out_transcript_ID(meta)

    known = 0
    ISM = 0
    prefix_ISM = 0
    suffix_ISM = 0
    NIC = 0
    NNC = 0
    genomic = 0
    antisense = 0
    intergenic = 0

    if "ISM_transcript" in meta: ISM = 1
    if "ISM-prefix_transcript" in meta: prefix_ISM = 1
    if "ISM-suffix_transcript" in meta: suffix_ISM = 1
    if "NIC_transcript" in meta: NIC = 1
    if "NNC_transcript" in meta: NNC = 1
    if "genomic_transcript" in meta: genomic = 1
    if "antisense_transcript" in meta: antisense = 1
    if "intergenic_transcript" in meta: intergenic = 1

    if sum([ISM, prefix_ISM, suffix_ISM, NIC, NNC, 
            genomic, antisense, intergenic]) == 0:
        known = 1

    return [transcript_ID, known, ISM, prefix_ISM, suffix_ISM, NIC, NNC, 
            genomic, antisense, intergenic]

def main():
    options = getOptions()
    gtf = options.gtf
    outprefix = options.outprefix

    # Open output files (BED and metadata)
    bed_file = open(outprefix + ".bed", 'w')
    metadata = open(outprefix + "_novelty.csv", 'w')
    meta_header = "\t".join(["transcript_ID", "known", "ISM", "prefix_ISM", "suffix_ISM", 
                             "NIC", "NNC", "genomic", "antisense", "intergenic"])
    metadata.write(meta_header + "\n")
   
    # Iterate over GTF and process transcripts
    with open(gtf, 'r') as f:
        for line in f:
            line = line.strip()
            entry = line.split("\t")
            
            # Skip things that are not transcripts
            if entry[2] != "transcript":
                continue
         
            # Get BED entry
            bed = create_BED_entry(entry) 
            bed_file.write("\t".join([str(x) for x in bed]) + "\n")

            # Get metadata entry
            meta = create_metadata_entry(entry)
            metadata.write("\t".join([str(x) for x in meta]) + "\n")

    bed_file.close()
    metadata.close()


if __name__ == '__main__':
    main()
