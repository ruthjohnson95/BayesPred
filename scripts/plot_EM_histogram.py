#!/usr/bin/env python

#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
from optparse import OptionParser
import sys
import numpy as np
import os
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="sim")
    parser.add_option("--results_file", dest="results_file")
    parser.add_option("--outdir", dest="outdir")
    (options, args) = parser.parse_args()

    # parse arguments
    name = options.name
    results_file = options.results_file
    outdir = options.outdir
    if outdir == None:
        outdir=""

    # read in file
    df = pd.read_csv(results_file, sep=' ')

    # remove zero entries 
    mu = df['mu']
    p = df['p']
    nonzero_inds=np.nonzero(p)

    mu = mu[nonzero_inds[0]]
    p = p[nonzero_inds[0]]

    df = {'mu': mu, 'p': p}
    df_nonzero = pd.DataFrame(data=df)


    # plot
    sns.set_style('whitegrid')

    df_nonzero.plot.bar(x='mu', y='p', legend=False)

    # save plot
    outfile=os.path.join(outdir, name + '.barplot.pdf')
    plt.title("Histogram of binned effect sizes - %s" % name, fontsize=15)
#    plt.xticks([])
#    plt.xlabel("0")

    plt.savefig(outfile)

if __name__== "__main__":
  main()
