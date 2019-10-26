# For each listed dataset, plot the read length distribution in a ridge plot

import sqlite3
import pandas as pd
from optparse import OptionParser
import os
import sys
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import joypy

def getOptions():
    parser = OptionParser()

    parser.add_option("--db", dest = "database",
        help = "TALON database", metavar = "FILE", type = "string")
    parser.add_option("--datasets", dest = "datasets",
        help = "Comma-delimited list of datasets to include", metavar = "FILE", type = "string")
    parser.add_option("--map", dest = "name_mapping",
        help = "Optional: a csv file mapping dataset IDsto other descriptors", metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outdir", help = "Directory for output files",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def get_read_lengths(database, datasets):
    """ Fetch the lengths of all reads from the specified datasets (list) and format
        in a pandas data frame """

    conn = sqlite3.connect(database)

    query = "SELECT dataset, read_length FROM observed WHERE dataset IN " + \
            format_for_IN(datasets)

    df = pd.read_sql_query(query, conn)
    conn.close()

    return df

def format_for_IN(l):
    """ Converts input to string that can be used for IN database query """

    if type(l) is tuple:
        l = list(l)
    if type(l) is str:
        l = [l]

    return "(" + ','.join(['"' + str(x) + '"' for x in l]) + ")"

def density_plot(data, outdir):
    """ Plot read length distribution"""

    fname = outdir + "/read_lengths.png"

    g = sns.FacetGrid(data, col='dataset')
    g.map(sns.distplot, "read_length")

    g.savefig(fname, dpi = 600)

def density_plot_with_mapping(data, outdir):
    """ Plot read length distribution for each dataset"""

    fname = outdir + "/read_lengths.png"

    g = sns.FacetGrid(data, row = 'name', hue = 'celltype')
    g.map(sns.distplot, "read_length")

    g.savefig(fname, dpi = 600)

def plot_indiv_dist(data, name, xmax, curr_color, outdir):
    name = "-".join(name.split())
    fname = outdir + "/" + name + "_read_lengths.png"
    xlabel = "Read length (bp)"
    ylabel = "Density"
    N50_label = "N50 = %d" % N50(data["read_length"])

    style = dict(size=18, color='black')
    sns.set_context("paper", font_scale=1.5)
    ax = sns.distplot(data["read_length"], color = curr_color)
    ax.set(xlabel=xlabel, ylabel=ylabel)   
    ax.set(xlim=(0, xmax))
    ax.set(ylim=(0, 0.0016))
    ymin, ymax = ax.get_ylim()
    ax.text(xmax*1/2, ymax*6/8, N50_label, **style)
 
    plot = ax.get_figure()
    plot.savefig(fname, dpi = 600) 
    plt.close(plot)

def N50(lengths):
    ## Sort reads longest>shortest
    all_len=sorted(lengths, reverse=True)
    csum=np.cumsum(all_len)

    n2=int(sum(lengths)/2)

    # get index for cumsum >= N/2
    csumn2=min(csum[csum >= n2])
    ind=np.where(csum == csumn2)

    n50 = all_len[ind[0][0]]
    return n50

def main():
    options = getOptions()
    database = options.database
    datasets = (options.datasets).split(",")
    mapping_file = options.name_mapping
    outdir = options.outdir

    # First, get data frame of read lengths
    data = get_read_lengths(database, datasets) 

    if mapping_file != None:
        # Parse name mapping file
        mapping = pd.read_csv(mapping_file)
        data = data.merge(mapping, on=['dataset', 'dataset'])
        xmax = data['read_length'].max()
        for dataset in datasets:
            print(dataset)
            curr = data.loc[data['dataset'] == dataset]
            name = curr.iloc[0]['name']
            color = curr.iloc[0]['color']
            plot_indiv_dist(curr, name, xmax, color, outdir)   
     
    else:
        # Plot the results in a ridge plot
        density_plot(data, outdir)
    
    # Print summary statistics
    grouped = data.groupby('dataset')
    print(grouped['read_length'].agg([N50]))
    print(grouped.size())

if __name__ == '__main__':
    main()
