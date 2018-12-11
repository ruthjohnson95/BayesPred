#!/usr/bin/env python

from optparse import OptionParser
import logging
import sys
import numpy as np
import scipy.stats as st
import os
import pandas as pd
import matplotlib.pyplot as plt

BOUND_A = -.20
BOUND_B = .20

def area_under(beta_hat):
    x_range = np.linspace(BOUND_A, BOUND_B)
    y = np.empty(len(x_range))
    M = len(beta_hat)

    for i,x in enumerate(x_range):
        y[i] = np.sum(beta_hat < x)/float(M)

    return y


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="sim")
    parser.add_option("--gwas_file", dest="gwas_file", default="/Users/ruthiejohnson/Development/unity_mixture/data/test_identity.2018.txt")
    parser.add_option("--results_file", dest="results_file", default="/Users/ruthiejohnson/Development/unity_mixture/data/test_identity.2018.results")
    parser.add_option("--outdir", dest="outdir", default="/Users/ruthiejohnson/Development/mixture_unity")
    parser.add_option("--title", dest="title", default="test_title")

    (options, args) = parser.parse_args()

    name = options.name
    gwas_file = options.gwas_file
    results_file = options.results_file
    outdir = options.outdir
    title = options.title

    # read in beta hats and true betas
    df_gwas = pd.read_csv(gwas_file, sep=' ')
    beta_hats = np.asarray(df_gwas['BETA_STD'])
    beta_true = np.asarray(df_gwas['BETA_TRUE'])

    # read in estimated p-vec
    df_results = pd.read_csv(results_file, sep=' ')
    mu_vec = np.asarray(df_results['mu'])
    p_est = np.asarray(df_results['p'])

    y_beta_hat = area_under(beta_hats)
    y_beta_true = area_under(beta_true)

    x = np.linspace(BOUND_A, BOUND_B)

    # plot estimate of p
    cumulative_p_est = np.cumsum(p_est)
    plt.plot(mu_vec, cumulative_p_est, color='red', label='estimated dist')

    # plot curves
    plt.plot(x, y_beta_hat, color='green', label='GWAS dist')
    plt.plot(x, y_beta_true, color='blue', label='true beta dist')

    # label plot
    plt.title("Estimated effect size dist - %s" % title)
    plt.xlabel("Effect size")
    plt.ylabel("Prop of density")
    plt.legend()


    plt.show()

if __name__== "__main__":
  main()
