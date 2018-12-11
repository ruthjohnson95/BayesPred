#!/usr/bin/env python

from optparse import OptionParser
import logging
import sys
import numpy as np
import scipy.stats as st
import os
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = OptionParser()
    parser.add_option("--gwas_file", dest="gwas_file", default="/Users/ruthiejohnson/Development/unity_mixture/data/test_identity.2018.txt")
    parser.add_option("--window", dest="window", default=10)
    (options, args) = parser.parse_args()

    gwas_file = options.gwas_file
    window = int(options.window)

    # read in beta hats and true betas
    df_gwas = pd.read_csv(gwas_file, sep=' ')
    beta_true = np.asarray(df_gwas['BETA_TRUE'])
    beta_hat = np.asarray(df_gwas['BETA_STD'])

    # prune by taking
    beta_true_prune = beta_true[0::window]
    beta_hat_prune = beta_hat[0::window]

    df = {'BETA_TRUE': beta_true_prune, 'BETA_STD': beta_hat_prune}
    df_prune = pd.DataFrame(data=df)
    df_prune.to_csv(gwas_file, index=False, sep=' ')

if __name__== "__main__":
  main()
