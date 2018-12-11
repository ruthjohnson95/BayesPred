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
    parser.add_option("--gwas_file", dest="gwas_file")
    parser.add_option("--outdir", dest="outdir")
    (options, args) = parser.parse_args()

    # parse arguments
    name = options.name
    gwas_file = options.gwas_file
    outdir = options.outdir
    if outdir == None:
        outdir=""

    # read in file
    df = pd.read_csv(gwas_file, sep=' ')
    beta_tilde = np.asarray(df['BETA_STD'])

    # plot
    sns.set_style('whitegrid')
#    a = -.50
#    b = .50
#    sns.distplot(beta_tilde, hist=True, kde=True, rug=False, hist_kws={"range": [a,b]})
    sns.distplot(beta_tilde, hist=True, kde=True, rug=False)

    # save plot
    outfile=os.path.join(outdir, name + '.hist.pdf')
    plt.title("Histogram of GWAS effect sizes - %s" % name, fontsize=15)

    plt.savefig(outfile)

if __name__== "__main__":
  main()
