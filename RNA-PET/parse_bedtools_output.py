# Consolidates outputs from running Bedtools Intersect on transcript starts/ends
# against RNA-PET start and end sites.
# Output format: transcript_ID, RNA-PET support (0/1)

from optparse import OptionParser

def getOptions():
    parser = OptionParser()

    parser.add_option("--starts", dest = "starts",
        help = "Bedtools intersect output for the starts", metavar = "FILE", type = "string")
    parser.add_option("--ends", dest = "ends",
        help = "Bedtools intersect output for the ends", metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def make_transcript_PET_dict(infile):
    """ Given a bedtools intersect file, this function creates a dictionary
        mapping each transcript_ID to a set containing the RNA-PET IDs that
        it matched to. """

    transcript_2_pet = {}   

    with open(infile, 'r') as f:
        for line in f:
            line = line.strip()
            entry = line.split("\t")

            transcript_ID = entry[3]
            rna_pet_ID = entry[9]

            if transcript_ID in transcript_2_pet:
                transcript_2_pet[transcript_ID].add(rna_pet_ID)
            elif rna_pet_ID != "-1":
                transcript_2_pet[transcript_ID] = set()
                transcript_2_pet[transcript_ID].add(rna_pet_ID)
            else:
               transcript_2_pet[transcript_ID] = set()
            
    return transcript_2_pet


def main():
    options = getOptions()
    start_file = options.starts
    end_file = options.ends
    outprefix = options.outprefix

    # Create dictionaries containing one key per input transcript pointing to
    # the set of all of its start or end RNA-PET matches (empty if none)
    start_matches = make_transcript_PET_dict(start_file)
    end_matches = make_transcript_PET_dict(end_file)

    # Now iterate over all of the transcripts and determine whether each one
    # had a successful start-end match pair or not
    o = open(outprefix + "_RNA-PET_results.csv", 'w')
    o.write(",".join(["transcript_ID", "RNA_PET_support"]) + "\n")
    transcript_IDs = list(start_matches.keys())
    for transcript in transcript_IDs:
        curr_starts = start_matches[transcript]
        curr_ends = end_matches[transcript]

        # If the union of the sets is smaller than the sum of their lengths,
        # that means that there was at least one start-end match
        if len(curr_starts | curr_ends) < len(curr_starts) + len(curr_ends):
            pet_support = "1"
        else:
            pet_support = "0"
        
        o.write(",".join([transcript, pet_support]) + "\n")

    o.close()
  

if __name__ == '__main__':
    main()
