# Convert FLAIR matrix to TALON-style abundance file
# usage: format_flair_matrix_like_talon.py <infile> <outfile> 


import sys

infile = sys.argv[1]
outfile = sys.argv[2]


o = open(outfile, 'w')
colnames = ["gene_ID", "transcript_ID", "annot_gene_id", 
            "annot_transcript_id", "annot_gene_name", "annot_transcript_name",
            "n_exons", "length", "gene_novelty", "transcript_novelty", 
            "ISM_subtype"]

line_num = 0
with open(infile, 'r') as f:
    for line in f:
        line = line.strip().split()

        if line_num == 0:
            sample_names = line[1:]
            o.write("\t".join(colnames + sample_names) + "\n")
        else:
            # Parse out IDs
            ids = line[0].split("_")
            gene_ID = ids[-1]
            transcript_ID = "_".join(ids[0:-1])

            entry = {k:"NA" for k in colnames}
            entry["gene_ID"] = entry["annot_gene_ID"] = entry["annot_gene_name"] = gene_ID
            entry["transcript_ID"] = entry["annot_transcript_ID"] = entry["annot_transcript_name"] = transcript_ID
            entry["counts"] = line[1:]

            # Determine novelty
            if "ENSG" in gene_ID:
                entry["gene_novelty"] = "Known"
            else:
                entry["gene_novelty"] = "Novel"

            if "ENST" in transcript_ID:
                entry["transcript_novelty"] = "Known"
            else:
                entry["transcript_novelty"] = "Novel"
             
            # Write to file
            outstr = "\t".join([entry["gene_ID"], entry["transcript_ID"],
                                entry["annot_gene_ID"], entry["annot_transcript_ID"],
                                entry["annot_gene_name"], entry["annot_transcript_name"],
                                entry["n_exons"], entry["length"], 
                                entry["gene_novelty"], entry["transcript_novelty"],
                                entry["ISM_subtype"]] + entry["counts"]) + "\n"
            o.write(outstr)
        line_num += 1
        

o.close()
