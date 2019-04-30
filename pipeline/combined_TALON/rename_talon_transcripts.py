# This script can be used to remap the TALON-issued gene and transcript names
# to a different prefix than the one used to initilize the TALON database.

import sqlite3
from sqlite3 import Error
from optparse import OptionParser
import os

def getOptions():
    parser = OptionParser()
    parser.add_option("--db", dest = "database",
        help = "TALON database",
        metavar = "FILE", type = "string")
    parser.add_option("--p", dest = "prefix",
        help = "New prefix for novel genes/transcripts",
        metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outname",
        help = "Outname for updated database (should end in '.db')",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def change_name_prefix(cursor, name):
    """ Go into the run_info section of the database and change the name prefix
    """
    update_query = 'UPDATE "run_info" SET "value" = ? WHERE "item" = ?'
    cursor.execute(update_query, [name, "idprefix"])

    return

def get_novel_gene_ids(cursor):
    """ Fetch IDs of all novel genes in the database """
    query = """SELECT DISTINCT(gene_ID) FROM observed
                   LEFT JOIN gene_annotations AS ga ON ga.ID = observed.gene_ID
                   WHERE (ga.attribute = 'gene_status' AND ga.value = 'NOVEL')
            """
    cursor.execute(query)
    novel_genes = [x[0] for x in cursor.fetchall()]
    return novel_genes

def get_novel_transcript_ids(cursor):
    """ Fetch IDs of all novel transcripts in the database """
    query = """SELECT DISTINCT(transcript_ID) FROM observed
                   LEFT JOIN transcript_annotations AS ta ON ta.ID = observed.transcript_ID
                   WHERE (ta.attribute = 'transcript_status' AND ta.value = 'NOVEL')
            """
    cursor.execute(query)
    novel_transcripts = [x[0] for x in cursor.fetchall()]
    return novel_transcripts

def update_gene_names(cursor, genes, prefix, n_places):
    """ Update TALON-issued gene names """

    for gene in genes:
        # Create the name
        gene_ID_str = str(gene).zfill(n_places)
        gene_name = prefix + "G" + gene_ID_str

        # Update the database
        update_query = """UPDATE "gene_annotations" SET "value" = ? 
                             WHERE "ID" = ? 
                             AND "attribute" = ? 
                             AND "source" = ? """
        cursor.execute(update_query, [gene_name, gene, "gene_name", "TALON"])
        cursor.execute(update_query, [gene_name, gene, "gene_id", "TALON"])

    return

def update_transcript_names(cursor, transcripts, prefix, n_places):
    """ Update TALON-issued transcript names """
    for transcript in transcripts:
        # Create the name
        ID_str = str(transcript).zfill(n_places)
        transcript_name = prefix + "T" + ID_str

        # Update the database
        update_query = """UPDATE "transcript_annotations" SET "value" = ?
                             WHERE "ID" = ?
                             AND "attribute" = ?
                             AND "source" = ? """
        cursor.execute(update_query, [transcript_name, transcript, "transcript_name", "TALON"])
        cursor.execute(update_query, [transcript_name, transcript, "transcript_id", "TALON"])

    return

def main():

    options = getOptions()
    orig_db = options.database
    new_db = options.outname
    prefix = options.prefix
    n_places = 9

    # Copy the database to the new name 
    os.system("cp %s %s" % (orig_db, new_db))

    # Connect to new database
    conn = sqlite3.connect(new_db)
    cursor = conn.cursor()
    conn.row_factory = sqlite3.Row

    # Update name in table
    change_name_prefix(cursor, prefix) 

    # Perform a query to grab the ID of every novel gene
    genes = get_novel_gene_ids(cursor)

    # Perform a query to grab the ID of every novel transcript
    transcripts = get_novel_transcript_ids(cursor)    

    # Now iterate over the genes and transcripts to update the names
    update_gene_names(cursor, genes, prefix, n_places) 
    update_transcript_names(cursor, transcripts, prefix, n_places)

    conn.commit()
    conn.close()

if __name__ == '__main__':
    main()
