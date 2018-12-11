#!/usr/bin/env python

from optparse import OptionParser
import logging
import sys
import numpy as np
import scipy.stats as st
import os
import pandas as pd
from sklearn.metrics import r2_score as r2_score

def sim_genotype_matrix(p, M, N):
    x = np.empty((N, M))
    for i in range(0, N):
        # x is individuals by snps (NxM)
        x_i = st.binom.rvs(2, p, size=M)
        x[i, :] = x_i
    return x

def calc_prs(x, weights):
    N = x.shape[0]
    prs = np.empty((N,1))

    for i in range(0, N):
        prs[i] = np.sum(np.multiply(x[i,:], weights))

    return prs


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="sim")
    parser.add_option("--M", dest="M")
    parser.add_option("--N", dest="N")
    parser.add_option("--maf", dest="maf")
    parser.add_option("--gwas_file", dest="gwas_file")
    parser.add_option("--outdir", dest="outdir", default="/Users/ruthiejohnson/Development/mixture_unity")

    (options, args) = parser.parse_args()

    M = int(options.M)
    N = int(options.N)
    maf = float(options.maf)
    outdir = options.outdir
    gwas_file = options.gwas_file

    # simulate genotypes
    x = sim_genotype_matrix(maf, M, N)

    # get weights from file
    df = pd.read_csv(gwas_file, sep=' ')
    true_weights = np.asarray(df['BETA_TRUE'])
    est_weights = np.asarray(df['WEIGHTS'])

    # compute weighted sum as prs
    true_prs = calc_prs(x, true_weights)
    est_prs = calc_prs(x, est_weights)

    # use r^2 as accuracy measurement
    r2 = r2_score(true_prs, est_prs)

    print("PRS r2: %.4f" % r2)


if __name__== "__main__":
  main()
