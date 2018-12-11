#!/usr/bin/env python

from optparse import OptionParser
import logging
import sys
import numpy as np
import scipy.stats as st
import os
import pandas as pd

# global variables
#logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)


def print_header(outfile, name, gwas_file, bins, ldsc_h2, N, outdir):

  header="BayesM (no LD)"
  print_func("%s" % header, outfile)
  print_func("Name: %s" % name, outfile)
  print_func("GWAS file: %s" % gwas_file, outfile)
  print_func("Bins: %d" % bins, outfile)
  print_func("LDSC h2: %.4f" % ldsc_h2, outfile)
  print_func("Sample size: %d" % N, outfile)
  print_func("Outdir: %s" % outdir, outfile)

  return

# prints both to console and to outfile with file descriptor f
def print_func(line, f):
    print(line)
    sys.stdout.flush()
    f.write(line)
    f.write('\n')
    return

def log_likelihood(beta_tilde, gamma, sigma_e, W):
    M = len(beta_tilde)
    mu = np.multiply(W, gamma)
    cov = np.multiply(np.eye(M), sigma_e)
    st.multivariate_normal.logpdf(x=beta_tilde, mean=mu, cov=cov)
    return log_like


def initialize_p(K):
    # random draw to initialize values
    p_init = np.random.dirichlet([1]*K,1)
    p_init = p_init.ravel()
    #p_init = [.5, .5]
    return p_init

def EM(p_t, mu_vec, sigma_vec, beta_tilde, sigma_e, its, f):
    K = len(p_t)
    M = len(beta_tilde)

    # print initial value
    p_init_string=""
    for p in p_t:
        p_init_string+= (str(p)+' ')

    print_func("Intial value for p: %s" % p_init_string, f)

    for i in range(0, its):
        # Expectation step
        r_MK = np.empty((M,K))

        for m in range(0, M):
            z = np.empty(K)
            for k in range(0, K):
                like = st.norm.pdf(beta_tilde[m], mu_vec[k], sigma_vec[k] + sigma_e)

                if np.nan == like:
                    like = 0
                    print_func("Encountered Nan!", f)
                z[k] = p_t[k]*like

            for k in range(0, K):
                if np.sum(z) > 0:
                    r_mk = z[k]/np.sum(z)
                else:
                    r_mk = 0

                r_MK[m,k] = r_mk

        # Maximization step
        for k in range(0, K):
            if np.sum(r_MK) > 0:
                p_t[k] = np.sum(r_MK[:,k])/np.sum(r_MK)
            else:
                p_t[k] = 0

        p_t_string = ""
        for p in p_t:
            p_t_string+=(str(p)+' ')

        print_func("Iteration %d: %s" % (i, p_t_string), f)


    return p_t


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="sim")
    parser.add_option("--gwas_file", dest="gwas_file")
    parser.add_option("--mu_vec", dest="mu_vec",)
    parser.add_option("--sigma_vec", dest="sigma_vec")
    parser.add_option("--ldsc_h2", dest="ldsc_h2")
    parser.add_option("--bins", dest="bins")
    parser.add_option("--ld_half_file", dest="ld_half_file")
    parser.add_option("--N", dest="N", default=100000)
    parser.add_option("--seed", dest="seed", default=100)
    parser.add_option("--outdir", dest="outdir", default="/Users/ruthiejohnson/Development/mixture_unity")
    parser.add_option("--its", dest="its", default=500)
    (options, args) = parser.parse_args()

    # parse command line args
    seed = int(options.seed)
    its = int(options.its)
    N = int(options.N)
    name = options.name
    gwas_file = options.gwas_file
    outdir = options.outdir
    ldsc_h2 = float(options.ldsc_h2)

    if options.bins is None:
        mu_vec = [float(item) for item in options.mu_vec.split(',')]
        sigma_vec = [float(item) for item in options.sigma_vec.split(',')]
    else:
        bins = int(options.bins)

        # get LHS of bins
        a = -.20
        b = .20
        step = (b-a)/float(bins)
        sigma_k = ((step*.5)/float(3))**2
        sigma_vec = np.repeat(sigma_k, bins)
        mu_vec = np.empty(bins)
        mu_vec_string = ""
        for k in range(bins):
            if k == 0:
                mu_vec[k] = a + step*.50
            else:
                mu_vec[k] = mu_vec[k-1] + step
            mu_vec_string+= (str(mu_vec[k]) +' ')

    # print header
    outfile=os.path.join(outdir, name+'.'+str(seed)+'.bayesM.log')
    f = open(outfile, 'w')
    print_header(f, name, gwas_file, bins, ldsc_h2, N, outdir)

    print_func("Mean of bins: %s" % mu_vec_string, f)

    print_func("Var of bins: %.4g" % sigma_vec[0], f)

    print_func("Reading in gwas file: %s" % gwas_file, f)

    df = pd.read_csv(gwas_file, sep=' ')

    ld_half_file = options.ld_half_file

    if ld_half_file is not None:
        beta_tilde = np.asarray(df['BETA_STD_I'])
    else:
        print_func("Assuming no LD", f)
        beta_tilde = np.asarray(df['BETA_STD'])

    # set seed
    np.random.seed(seed)

    if ld_half_file is not None:
        W = np.loadtxt(ld_half_file)
        print_func("Using ld half file: %s" % ld_half_file, f)

    # get meta values from data
    K = len(mu_vec)
    M = beta_tilde.shape[0]

    # intialize values of chain
    print_func("initializing p", f)
    p_init = initialize_p(K)

    # calculate sigma_e
    sigma_e = (1-ldsc_h2)/N

    p_est = EM(p_init, mu_vec, sigma_vec, beta_tilde, sigma_e, its, f)

    p_est_str=""
    for p in p_est:
        p_est_str += (str(p) + ' ')

    print_func("Estimate of proportions: %s" % p_est_str, f)

    # put results in a dataframe
    df = {'mu': mu_vec, 'p': p_est}
    results_df = pd.DataFrame(data=df)
    results_file = os.path.join(outdir, name +'.'+str(seed)+'.results')
    results_df.to_csv(results_file, index=False, sep=' ')

    f.close()


if __name__== "__main__":
  main()
