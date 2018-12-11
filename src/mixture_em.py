#!/usr/bin/env python

from optparse import OptionParser
import logging
import sys
import numpy as np
import scipy.stats as st
import os
from scipy.special import logsumexp

# global variables
logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

VARY_MU = False
VARY_SIGMA = False
MAX_ITS = 10

def precompute(mu_vec, W):
    K = len(mu_vec)
    M = W.shape[0]

    nu = np.empty((M, K))

    for k, mu_k in enumerate(mu_vec):
        mu_k_vec = np.repeat(mu_k, M)
        nu[:,i] = np.matmul(W, mu_k_vec)

    return nu


def initialize(K):
    # random draw to initialize values
    p_init = np.random.dirichlet([1]*K,1)
    return p_init


def expectation_step(p_K, nu, sigma_K, beta_hat, N):
    M = nu.shape[0]
    K = nu.shape[1]

    r_MK = np.empty((M,K))

    for k in range(0, K):
        for m in range(0,M):
            M_k = np.multiply(M, p_K)
            sigma_g = np.sum(np.multiply(M_k, sigma_K))
            sigma_e = (1-sigma_g)/float(N)

            mu_mk = nu[m,k]
            sigma_mk = sigma_e + sigma_K[k]
            numerator = p_K[k] + st.norm.pdf(beta_hat[m], mu_mk, sigma_mk)
            denominator = 0
            for l in range(0, K):
                mu_ml = nu[m,l]
                sigma_ml = sigma_e + sigma_K[l]
                denominator += p_K[l] * st.norm.pdf(beta_hat[m], mu_ml, sigma_ml)

            r_MK[m,k] = numerator/denominator

    return r_MK


def maximization_step(r_MK):
    M = r_MK.shape[0]
    p_K = np.divide(np.sum(r_MK, axis=1), M)
    return p_K


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="sim")
    parser.add_option("--gwas_file", dest="gwas_file", default="sim.100.txt")
    parser.add_option("--mu_vec", dest="mu_vec", nargs='+', default="-.10 0 .10")
    parser.add_option("--sigma_vec", dest="sigma_vec", nargs="+", default=".001 .001 .001")
    parser.add_option("--bins", dest="bins")
    parser.add_option("--ld_half_file", dest="ld_half_file")
    parser.add_option("--N", dest="N", default=100000)
    parser.add_option("--seed", dest="seed", default=100)
    parser.add_option("--outdir", dest="outdir", default="/Users/ruthiejohnson/Development/mixture_unity")
    parser.add_option("--precompute", dest="precompute", default='y')
    (options, args) = parser.parse_args()

    """
    if mu_vec is None and bins is None:
        global VARY_MU
        VARY_MU = True
        logging.info("Estimating mixture means: mu_0,...,mu_K")

    if sigma_vec is None and bins is None:
        global VARY_SIGMA
        VARY_SIGMA = True
        logging.info("Estimating mixture variances: sigma_0,...,sigma_K")
    """

    logging.info("Estimating mixing proportions: p_0,...,p_K")

    ld_half_file = options.ld_half_file
    W = np.loadtxt(ld_half_file)
    logging.info("Using ld half file: %s" % ld_half_file)

    if precompute == 'y':
        logging.info("Pre-computing terms")
        nu = precompute(mu_vec, W)

    p_K = initialize(K)

    for i in range(0, MAX_ITS):

        # expectation step
        r_MK = expectation_step(p_K, nu, sigma_K, beta_hat, N)

        # maximization step
        p_K = maximization_step(r_MK)

    # end EM-loop

    print(p_K)

    return


if __name__== "__main__":
  main()
