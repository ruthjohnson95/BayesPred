#!/usr/bin/env python

from optparse import OptionParser
import logging
import sys
import numpy as np
import scipy.stats as st
import os
import pandas as pd

"""
    prints both to console and to outfile with file descriptor f
"""
def print_func(line, f):
    print(line)
    sys.stdout.flush()
    f.write(line)
    f.write('\n')
    return

"""
    Computes log-likelihood of GWAS effect sizes given latent variables
"""
def log_likelihood(beta_tilde, gamma, sigma_e, W):
    M = len(beta_tilde)
    mu = np.multiply(W, gamma)[0]
    cov = np.multiply(np.eye(M), sigma_e)
    log_like = st.multivariate_normal.logpdf(x=beta_tilde, mean=mu, cov=cov)
    return log_like


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="simulated")
    parser.add_option("--ld_half_file", dest="ld_half_file")
    parser.add_option("--gwas_file", dest="gwas_file")
    parser.add_option("--ldsc_h2", dest="ldsc_h2")
    parser.add_option("--N", dest="N")
    parser.add_option("--seed", dest="seed", default=100)
    parser.add_option("--outdir", dest="outdir", default="/Users/ruthiejohnson/Development/unity_mixture/data")

    (options, args) = parser.parse_args()

    logging.info("Computing true log-likelihood")

    # parse input options
    outname = options.name
    outdir = options.outdir
    ldsc_h2 = float(options.ldsc_h2)
    N = int(options.N)
    seed = int(options.seed)
    gwas_file = options.gwas_file

    # save to outfile
    outfile = os.path.join(outdir, outname + '.' + str(seed) + '.like.txt')
    f = open(outfile, 'w')

    # read in LD
    ld_half_file = options.ld_half_file
    W = np.loadtxt(ld_half_file)

    # read in true betas
    df = pd.read_csv(gwas_file, sep=' ')
    true_betas = np.asarray(df['BETA_TRUE'])
    beta_tilde = np.asarray(df['BETA_STD_I'])

    # calculate sigma_e
    sigma_e = (1-ldsc_h2)/float(N)
    log_like = log_likelihood(beta_tilde, true_betas, sigma_e, W)

    # print to file
    print_func("Log_like: %.4g" % log_like, f)

    f.close()


if __name__== "__main__":
  main()
