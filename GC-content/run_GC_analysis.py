
from optparse import OptionParser
import os
import sys
from Bio import SeqIO
from Bio.SeqUtils import GC
import gzip
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

def getOptions():
    parser = OptionParser()

    parser.add_option("--fasta", dest = "fasta",
        help = "Fasta file of GENCODE transcript sequences", metavar = "FILE", type = "string")
    parser.add_option("--genes", dest = "genes",
        help = "CSV file of gene names and whether PacBio detected them or not", metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outdir", help = "Directory for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def get_gene_lists(infile):
    """ Read in gene list file and group genes into two sets based on whether
        each was detected in PacBio or not """

    detected = set()
    underdetected = set()

    with open(infile, 'r') as f:
        for line in f:
            line = line.strip()

            gene,detection = line.split(",")
            if detection == "FALSE":
                underdetected.add(gene)
            elif detection == "TRUE":
                detected.add(gene)
            else:
                sys.exit("Unexpected value in gene list: '" + detection + "'")

    return detected,underdetected

def compute_all_GCs(fasta, detected_genes, underdetected_genes):
    """ For each fasta transcript:
          1) Extract gene name
          2) Check if gene name appears in detected/underdetected gene lists.
             If not, discard. If yes, proceed.
          3) Compute GC content of sequence
          4) Record gene name, transcript ID, detection status, and GC content
             in pandas table.
    """
    gene_names = []
    transcript_IDs = []
    detection_status = []
    GC_content = []

    try:
        with gzip.open(fasta, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                split_ID = (record.id).split("|")
                gene_name = split_ID[5]
                transcript_ID = split_ID[0]

                if gene_name not in detected_genes and \
                   gene_name not in underdetected_genes:
                    continue
                if gene_name in detected_genes:
                    detection_status.append("Yes")
                elif gene_name in underdetected_genes:
                    detection_status.append("No")

                gene_names.append(gene_name)
                transcript_IDs.append(transcript_ID)
                GC_content.append(GC(record.seq))
               
    except Exception as e:
        print(e)
        sys.exit("Problem reading fasta sequence file. Expecting gzipped file")

    # Convert lists into a pandas data frame
    df = pd.DataFrame({"gene_name": gene_names, 
                       "transcript_ID": transcript_IDs, 
                       "detected": detection_status, 
                       "GC": GC_content})
   
    return df 

def GC_violin_plot(data, outdir):
    """ Plot GC content for detected/underdetected genes """

    fname = outdir + "/GC_plot.png"
    xlabel = "Gene detected in PacBio"
    ylabel = "Transcript GC content"

    sns.set_context("paper", font_scale=1.5)
    ax = sns.violinplot(data = data, x="detected", y="GC", palette="deep")
    ax.set(xlabel=xlabel, ylabel=ylabel)

    plot = ax.get_figure()
    plot.savefig(fname, dpi = 600)

def main():
    options = getOptions()
    fasta = options.fasta
    genes = options.genes
    outdir = options.outdir

    # Create two gene lists: detected and underdetected
    detected,underdetected = get_gene_lists(genes)    

    # Read fasta file and compute GC content of entries 
    GC_table = compute_all_GCs(fasta, detected, underdetected)

    # Plot GC content of detected vs underdetected genes
    GC_violin_plot(GC_table, outdir)
    grouped = GC_table.groupby('detected')
    print(grouped['GC'].agg([np.sum, np.mean, np.std]))
    print(grouped.size())



if __name__ == '__main__':
    main()
