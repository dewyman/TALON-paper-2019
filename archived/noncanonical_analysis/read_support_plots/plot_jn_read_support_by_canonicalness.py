# Take a splice junction file and plot the amount of read support for
# both canonical and noncanonical splice junctions as a violin plot. 

from optparse import OptionParser
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
sns.set(style="whitegrid")

def getOptions():
    parser = OptionParser()
   
    parser.add_option("--sjs", dest = "sj_file",
        help = "File of splice junctions", metavar = "FILE", type = "string")
    parser.add_option("--ymax", dest = "ymax",
        help = "Max y value to show", metavar = "FILE", type = "string")
    parser.add_option("--o", dest = "outprefix", help = "Prefix for output file",
        metavar = "FILE", type = "string")

    (options, args) = parser.parse_args()
    return options

def label_canon(row):
    """ Return "Canonical" if row is a canonical junction, "Noncanonical"
        otherwise. """
    if row['motifcode'] == 0 :
      return 'Noncanonical'
    else:
      return 'Canonical'

def main():

    options = getOptions()
    sj_file = options.sj_file
    ymax = options.ymax
    fname = options.outprefix + "_read_support_by_jn_canonicalness.png"

    # Read data into pandas df
    junctions = pd.read_csv(sj_file, sep = '\t', header=None)
    junctions.columns = ['chrom', 'start', 'end', 'strand', 'motifcode',
                         'annot', 'nReads','mult','mx']

    # Create new Canonical/Noncanonical column
    junctions['canon'] = junctions.apply (lambda row: label_canon(row), axis=1)
    print(max((junctions.loc[junctions.canon == "Noncanonical"]).nReads))

    # Plot
    if ymax == None:
        ymax = max(junctions['nReads'])
        ymax = ymax + ymax*0.05
    else:
        ymax = int(ymax) 

    ax = sns.violinplot(x="canon", y="nReads", data=junctions, 
                        hue="canon", cut=0,
                        saturation = 1, inner = "box", linewidth = 1)
    ax.set(ylim=(0, ymax))
    ax.set(xlabel='Splice junction type', ylabel='Number of mapped long reads supporting junction')

    plt.savefig(fname, dpi = 600)
    plt.close('all')

    # Write stats to file
    df = junctions.groupby('canon')['nReads'].describe()[["mean"]]
    outfile = options.outprefix + "_read_support_by_jn_canonicalness.txt"
    df.to_csv(outfile)

if __name__ == '__main__':
    main()
