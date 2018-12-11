#!/usr/bin/env python

from optparse import OptionParser
import logging
import sys
import numpy as np
import scipy.stats as st
import os
import pandas as pd


def simulate_ld(M, coef):
    d = np.ones(M)*coef
    V = np.diag(d)

    # decreasing values along the diagonal
    for m in range(0, M):
        center = V[m,m]
        power = 0
        for l in range(m, M):
            V[m, l] = center ** power
            power += 1
    # reflect to bottom triangluar
    for i in range(M):
        for j in range(i, M):
            V[j][i] = V[i][j]

    # ensure pos-semi ef
    V[:] = truncate_matrix(V)

    return V


def truncate_eigenvalues(d):
    M = len(d)

    # order evaules in descending order
    d[::-1].sort()

    #running_sum = 0
    d_trun = np.zeros(M)

    # keep only positive evalues
    for i in range(0,M):
        if d[i] > 0:
            # keep evalue
            d_trun[i] = d[i]

    return d_trun


def truncate_matrix(V):
    # make V pos-semi-def
    d, Q = np.linalg.eigh(V, UPLO='U')

    # reorder eigenvectors from inc to dec
    idx = d.argsort()[::-1]
    Q[:] = Q[:, idx]

    # truncate small eigenvalues for stability
    d_trun = truncate_eigenvalues(d)

    # mult decomp back together to get final V_trunc
    M1 = np.matmul(Q, np.diag(d_trun))
    V_trun = np.matmul(M1, np.matrix.transpose(Q))

    return V_trun


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="simulated")
    parser.add_option("--M", dest="M", default=10)
    parser.add_option("--coef", dest="coef", default=0)
    parser.add_option("--seed", dest="seed", default=100)
    parser.add_option("--outdir", dest="outdir", default="/Users/ruthiejohnson/Development/unity_mixture/data")

    (options, args) = parser.parse_args()

    # parse input options
    outname = options.name
    M = int(options.M)
    seed = int(options.seed)
    np.random.seed(seed)
    outdir = options.outdir
    coef = float(options.coef)

    V = simulate_ld(M, coef)

    # save to outfile
    outfile = os.path.join(outdir, outname + '_' + str(coef) + '.' + str(M) + '.txt')
    np.savetxt(outfile, V)

    logging.info("DONE...simulated LD can be found at: %s" % outfile)


if __name__== "__main__":
  main()
