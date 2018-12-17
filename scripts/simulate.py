#!/usr/bin/env python

from optparse import OptionParser
import logging
import sys
import numpy as np
import scipy.stats as st
import os
import pandas as pd

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

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


def simulate_mixture(p_vec, mu_vec, sigma_vec, M, sigma_e, V):
    K = len(p_vec)
    k_matrix = np.random.multinomial(1, p_vec, size=M)

    # make large beta_k_matrix
    beta_k_matrix = np.zeros((M,K))

    # sample from K gaussians
    for k in range(0, K):
        beta_k_matrix[:,k] = st.norm.rvs(loc=mu_vec[k],scale=sigma_vec[k],size=M)

    # only pick betas from chosen components
    beta_k = np.multiply(beta_k_matrix,  k_matrix)

    print "True p"
    print np.divide(np.sum(k_matrix, axis=0), float(M))

    beta = np.sum(beta_k, axis=1)

    # sample GWAS effects
    mu = np.matmul(V, beta)
    cov = np.multiply(V, sigma_e)

    beta_hat = st.multivariate_normal.rvs(mu, cov)
    beta_true = beta

    df = {'BETA_STD': beta_hat, 'BETA_TRUE': beta_true, 'GAMMA_1': beta_k[:,0], 'GAMMA_2':beta_k[:,1], 'GAMMA_3':beta_k[:,2], 'GAMMA_4':beta_k[:,3]}
    beta_hat_df = pd.DataFrame(data=df)

    return beta_hat_df


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="sim")
    parser.add_option("--p_vec", dest="p_vec", default=".25 .25 .50")
    parser.add_option("--mu_vec", dest="mu_vec", default="-.10 0 .10")
    parser.add_option("--sigma_vec", dest="sigma_vec", default=".001 .001 .001")
    parser.add_option("--bins", dest="bins")
    parser.add_option("--ld_file", dest="ld_file")
    parser.add_option("--M", dest="M", default=10)
    parser.add_option("--N", dest="N", default=100000)
    parser.add_option("--sigma_g", dest="sigma_g", default=.50)
    parser.add_option("--seed", dest="seed", default=100)
    parser.add_option("--outdir", dest="outdir", default="/Users/ruthiejohnson/Development/mixture_unity")

    (options, args) = parser.parse_args()

    # parse input options
    if options.bins is None:
        p_vec= [float(item) for item in options.p_vec.split(',')]
        mu_vec = [float(item) for item in options.mu_vec.split(',')]
        sigma_vec = [float(item) for item in options.sigma_vec.split(',')]
    else:
        bins = int(options.bins)
        p_vec= [float(item) for item in options.p_vec.split(',')]

        if bins != len(p_vec):
            logging.info("Error: number of bins does not equal mixing proportions...exiting")
            exit(1)

        # get LHS of bins
        a = -.20
        b = .20
        step = (b-a)/float(bins)
        sigma_k = ((step*.5)/float(3))**2
        sigma_vec = np.repeat(sigma_k, bins)
        mu_vec = np.empty(bins)
        for k in range(bins):
            if k == 0:
                mu_vec[k] = a + step*.50
            else:
                mu_vec[k] = mu_vec[k-1] + step

        print "Sigma vec:"
        print sigma_vec
        print "mu vec:"
        print mu_vec


    outname = options.name
    ld_file = options.ld_file
    M = int(options.M)
    N = int(options.N)
    seed = int(options.seed)
    np.random.seed(seed)
    outdir = options.outdir
    sigma_g = float(options.sigma_g)

    # check that mixture components are same size
    if len(p_vec) == len(mu_vec) == len(sigma_vec):
        pass
    else:
        logging.info("ERROR: p, mu, sigma lists are NOT same size...exiting")
        exit(1)

    # simulate values
    print ld_file
    if ld_file is None:
        V = np.eye(M)
    else:
        try:
            V_raw = np.loadtxt(ld_file)
            # truncate matrix to make pos-semi def
            logging.info("Truncating matrix to ensure pos-semi-def")
            V = truncate_matrix(V_raw)

        except:
            logging.info("LD file does not exist...will simulate with no LD")
            V = np.eye(M)

    # calculate sigma_e
    sigma_e = (1-sigma_g)/N

    beta_hats = simulate_mixture(p_vec, mu_vec, sigma_vec, M, sigma_e, V)

    # save to outfile
    outfile = os.path.join(outdir, outname + '.' + str(seed) + '.txt')
    beta_hats.to_csv(outfile, index=False, sep=' ')

    logging.info("DONE...simulated files can be found at: %s" % outfile)


if __name__== "__main__":
  main()
