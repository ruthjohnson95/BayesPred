#/usr/bin/env python

"""
    This program models the effect size distribution of effect sizes from GWAS summary statistics
    through a mixture of Gaussians. The main parameter of interest is to estimate the proportion
    of SNPs that correspond to each mixing component.

    INPUT: standardized GWAS effect sizes
    OUTPUT: proportion of SNPs per kth bin (p_k)
"""

from optparse import OptionParser
import logging
import sys
import numpy as np
import scipy.stats as st
import os
import pandas as pd
import math
from scipy.misc import logsumexp as logsumexp
from sklearn.metrics import r2_score as r2_score
import cProfile, pstats, StringIO

# global variables
logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)
EXP_MAX = math.log(sys.float_info.max)
EXP_MIN = math.log(sys.float_info.min)
OPT = True

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


"""
    Computes denominator for posterior dist of mixture assignments
"""
def compute_km_vec(mu_km_vec, sigma_km_vec, mu_vec, sigma_vec, p_vec):
    sum = 0
    K = len(mu_km_vec)
    q_km_vec = np.empty(K)

    a_vec = np.empty(K)
    var_vec = np.empty(K)

    for k in range(0,K):
        sigma_k = sigma_vec[k]
        mu_k = mu_vec[k]
        var_term = np.sqrt(sigma_km_vec[k]/sigma_vec[k])
        temp_term_1 = .5 * 1 / float(sigma_km_vec[k]) * mu_km_vec[k] * mu_km_vec[k]
        temp_term_2 = .5* 1/ float(sigma_k) * mu_k * mu_k
        temp_term = temp_term_1 - temp_term_2

        a_vec[k] = temp_term[0]
        var_vec[k] = var_term[0]

        """
        #print temp_term
        if temp_term > EXP_MAX:
            print "EXP_MAX"
            #print temp_term
            #print "sigma_km: %.4g" % sigma_km_vec[k]
            #print "sigma_k: %.4g" % sigma_k
            temp_term = EXP_MAX
        elif temp_term < EXP_MIN:
            print "EXP_MIN"
            #print temp_term
            temp_term = EXP_MIN

        mean_term = np.exp(temp_term)


        q_km_vec[k] = p_vec[k]*var_term*mean_term
        """

    # calculate q_km in log form
    log_terms = np.empty(K)
    for k in range(0, K):
        if p_vec[k] == 0 or var_vec[k] == 0:
            print "log error!"
            exit(1)

        log_terms[k] = np.log(p_vec[k]) + np.log(var_vec[k]) + a_vec[k]

    q_log_denom = logsumexp(log_terms)

    for k in range(0, K):
        q_log_num = np.log(p_vec[k]) + np.log(var_vec[k]) + a_vec[k]
        q_km_vec[k] = np.exp(q_log_num - q_log_denom)

    return q_km_vec


"""
    Computes parameters for posterior distriubtion of mixture assignments (C_k)
"""
def compute_q_km(k, m, p_vec, mu_vec, sigma_vec, psi_m, A, gamma_t, C_t, sigma_e, W, beta_tilde):

    # holds means and variances across K components (NOT SNPs)
    K = len(mu_vec)
    mu_km_vec = np.empty((K,1))
    sigma_km_vec = np.empty((K,1))
    q_vec=np.zeros(K)

    for k in range(0, K):
        mu_k = mu_vec[k]
        sigma_k = sigma_vec[k]
        mu_km_vec[k], sigma_km_vec[k] = compute_mu_sigma_km_opt(m, k, mu_k, sigma_k, psi_m, A, gamma_t, C_t, sigma_e, W, beta_tilde)

    q_km_vec = compute_km_vec(mu_km_vec, sigma_km_vec, mu_vec, sigma_vec, p_vec)

    q_km_denom = np.sum(q_km_vec)
    for k in range(0, K):
        q_km_num = q_km_vec[k]
        q_km = q_km_num/q_km_denom
        q_vec[k] = q_km

    return q_vec


"""
    Computes sum for effient update
"""
def dp_term(m, A, gamma_k):
    sum = 0
    nonzero_inds = np.nonzero(gamma_k)[0]
    for i in nonzero_inds:
        if i != m:
            sum += A[i,m] * gamma_k[i]
    return sum


"""
    Computes mean and standard deviation for posterior distribution
    of effect sizes for kth component
"""
def compute_mu_sigma_km(m, k, mu_k, sigma_k, psi_m, A, gamma_t, C_t, sigma_e, W, beta_tilde):

    sigma_km = 1/(1/sigma_k + 1/sigma_e)

    # slow way
    W_m = W[:, m]
    gamma_C_t = np.sum(np.multiply(C_t, gamma_t), axis=1)

    r_m_1 = np.subtract(beta_tilde, np.matmul(W, gamma_C_t))
    r_m_2 = np.multiply(W_m, gamma_t[m, k])
    r_m = np.add(r_m_1, r_m_2)
    r_m_T = np.transpose(r_m)
    mu_km = sigma_km*((mu_k/sigma_k) + ((np.matmul(r_m_T, W_m))/sigma_e))

    # faster way
    #term1 = (sigma_km*mu_k)/sigma_k
    #term2 = sigma_km*sigma_e
    #dp = dp_term(m, A, gamma_C_t)
    #term3 = psi_m - dp
    #mu_km = term1 + term2*term3

    return mu_km, sigma_km


def compute_mu_sigma_km_opt(m, k, mu_k, sigma_k, psi_m, A, gamma_t, C_t, sigma_e, W, beta_tilde):

    sigma_km = 1/(1/sigma_k + 1/sigma_e)

    M = len(beta_tilde)

    W_m = W[:, m]

    # zero out last row
    #C_t[:,3] = np.zeros(M)

    #gamma_C_t = np.sum(np.multiply(C_t, gamma_t), axis=1)
    gamma_C_t = np.sum(np.multiply(C_t[:,:3], gamma_t[:,:3]), axis=1)

    nonzero_inds = np.nonzero(gamma_C_t)[0]

    sum = 0
    
    for i in nonzero_inds:
    #for i in range(0, M):
        a_im = A[i,m]
        if i!=m:
            # find nonzero element of C_t
            #ind = np.argmax(C_t[i,:])
            #sum += a_im * gamma_t[i, ind]
            sum += a_im * gamma_C_t[i]

    dp_term = psi_m - sum

    mu_km = sigma_km*((mu_k/sigma_k) + ((dp_term)/sigma_e))

    return mu_km, sigma_km


"""
    Main Gibbs sampler function
"""
def gibbs(p_init, gamma_init, C_init, mu_vec, sigma_vec, sigma_e, W, A, psi, beta_tilde, N, its, f):

    # get metadata
    K = len(mu_vec)
    M = len(beta_tilde)

    # make lists to hold samples
    p_list = []
    C_list = np.empty((M, K))

    weights = np.zeros(M)

    # start chain
    p_t = p_init
    C_t = C_init

    gamma_t = gamma_init

    # start sampler
    logging.info("Starting sampler")
    for i in range(0, its): # run for iterations
        for m in range(0, M): # loop through M SNPs

            # temp vector to hold effect size estimates
            gamma_temp = np.empty(K)

            # sample gamma_t
            for k in range(0, K): # loop through K mixture components
                # if SNP belongs to mixture K, otherwise effect is 0
                if C_t[m,k] == 1:
                    # compute posterior mean and variance
                    mu_km, sigma_km = compute_mu_sigma_km_opt(m, k, mu_vec[k], sigma_vec[k], psi[m], A, gamma_t, C_t, sigma_e, W, beta_tilde)

                    # sample effect sizes from the posterior
                    gamma_temp[k] = st.norm.rvs(mu_km, sigma_km)
                    #gamma_t[m,k] = st.norm.rvs(mu_km, sigma_km)
                else:
                    #gamma_t[m,k] = 0
                    gamma_temp[k] = 0

            # only update gamma_t after all K mixtures have been sampled
            gamma_t[m,:] = gamma_temp

            # sample mixture assignments
            q_km = compute_q_km(k, m, p_t, mu_vec, sigma_vec, psi[m], A, gamma_t, C_t, sigma_e, W, beta_tilde)
            C  = st.multinomial.rvs(n=1, p=q_km, size=1)
            C = C.ravel()
            C_t[m,:] = C

            # end loop through K clusters
        # end loop through SNPs


            # end loop through K clusters
        # end loop through SNPs

        alpha = np.add(np.sum(C_t, axis=0), np.ones(K))

        p_t = st.dirichlet.rvs(alpha)
        p_t = p_t.ravel()
        p_list.append(p_t)

        gamma_C_t = np.sum(np.multiply(C_t, gamma_t),axis=1)

        # speedup
        #log_like = log_likelihood(beta_tilde, gamma_C_t, sigma_e, W)
        log_like = 0

        p_t_string = ""
        for p in p_t:
            p_t_string+=(str(p)+' ')

        if i%50 == 0:
            print_func("Iteration %d: %s" % (i, p_t_string), f)
            #print_func("Iteration %d (log-like): %4g" % (i, log_like), f)


        # compute weight for iteration
        BURN = its/4
        if i >= BURN:
            weights[:] = np.add(weights, gamma_C_t)

    # end loop

    p_est = np.mean(p_list[BURN:], axis=0)

    C_est = np.divide(np.sum(C_t, axis=0), float(M))

    return p_est, C_est, weights


"""
    Pre-computes transformed effect sizes and matrix inverse
    for efficient update
"""
def precompute_terms(W, beta_tilde, name, outdir):
    A = np.matmul(np.transpose(W), W)
    M = W.shape[0]
    psi = np.empty((M,1))
    for m in range(0,M):
        psi[m] = np.matmul(np.transpose(beta_tilde), W[:,m])

    # save for future use
    A_file = os.path.join(outdir, name+'.A')
    psi_file = os.path.join(outdir, name+'.psi')
    np.save(A_file, A)
    np.save(psi_file, psi)
    return A, psi


"""
    Initializes mixture proportions
"""
def initialize_p(K):
    # random draw to initialize values
    p_init = np.random.dirichlet([10]*K,1)
    p_init = p_init.ravel()

    return p_init


"""
    Initialize effect sizes and mixture assignements for each SNP
"""
def initialize_C_gamma(p_init, mu_vec, sigma_vec, M):
    # create empty array to hold values
    K = len(mu_vec)
    gamma_init = np.empty((M,K))

    C_init = np.random.multinomial(n=1, pvals=p_init, size=M)

    for k in range(0, K):
        gamma_init[:, k] = st.norm.rvs(mu_vec[k], sigma_vec[k], size=M)

    return C_init, gamma_init


def main():
    parser = OptionParser()
    parser.add_option("--name", dest="name", default="sim")
    parser.add_option("--gwas_file", dest="gwas_file")
    parser.add_option("--mu_vec", dest="mu_vec",)
    parser.add_option("--sigma_vec", dest="sigma_vec")
    parser.add_option("--bins", dest="bins")
    parser.add_option("--ld_half_file", dest="ld_half_file")
    parser.add_option("--N", dest="N", default=100000)
    parser.add_option("--seed", dest="seed", default=100)
    parser.add_option("--outdir", dest="outdir", default="/Users/ruthiejohnson/Development/mixture_unity")
    parser.add_option("--precompute", dest="precompute", default='y')
    parser.add_option("--ldsc_h2", dest="ldsc_h2", default=0.50)
    parser.add_option("--its", dest="its", default=500)
    parser.add_option("--opt", dest="opt", default="y")
    (options, args) = parser.parse_args()

    # parse command line args
    seed = int(options.seed)
    its = int(options.its)
    N = int(options.N)
    name = options.name
    precompute = options.precompute
    gwas_file = options.gwas_file
    outdir = options.outdir
    ldsc_h2 = float(options.ldsc_h2)

    # optimze speedup
    opt = options.opt
    if opt != 'y':
        print "NOT using speedup"
        global OPT
        OPT = False
    else:
        print "Using speedup"


    # log file
    outfile=os.path.join(outdir, name+'.'+str(seed)+'.BayesPred.log')
    f = open(outfile, 'w')

    if options.bins is None and options.mu_vec is not None and options.sigma_vec is not None: # user provides mean/variances for mixture components
        mu_vec = [float(item) for item in options.mu_vec.split(',')]
        sigma_vec = [float(item) for item in options.sigma_vec.split(',')]
    elif options.bins is not None: # create evenly spaced bins and use mean/variance computed from each bin size
        bins = int(options.bins)
        a = -.10 # LHS of interval
        b = .10 # RHS of interval
        step = (b-a)/float(bins)
        sigma_k = ((step*.5)/float(3))**2
        sigma_vec = np.repeat(sigma_k, bins)
        mu_vec = np.empty(bins)
        mu_vec_string = ""

        # put mu vec into string for printing
        for k in range(bins):
            if k == 0:
                mu_vec[k] = a + step*.50
            else:
                mu_vec[k] = mu_vec[k-1] + step
            mu_vec_string+= (str(mu_vec[k]) +' ')
        print "mean vec: %s" % mu_vec_string
        print "sigma vec: "
        print sigma_vec

    else:
        logging.info("ERROR: user needs to specify mu/sigma or bins...exiting")
        exit(1)

    # read in transformed GWAS effect sizes
    logging.info("Reading in gwas file: %s" % gwas_file)
    df = pd.read_csv(gwas_file, sep=' ')
    beta_tilde = np.asarray(df['BETA_STD_I'])

    # set seed
    np.random.seed(seed)

    # read in LD
    ld_half_file = options.ld_half_file
    W = np.loadtxt(ld_half_file)
    logging.info("Using ld half file: %s" % ld_half_file)

    # precompute terms or load in
    if precompute == 'y':
        logging.info("Pre-computing terms")
        A, psi = precompute_terms(W, beta_tilde, name, outdir)
    else: # assumes files have already been made
        # check if files exist
        A_file = os.path.join(outdir, name+'.A.npy')
        psi_file = os.path.join(outdir, name+'.psi.npy')
        try:
            logging.info("Reading in pre-computed matrices")
            A = np.load(A_file)
            psi = np.load(psi_file)
        except:
            logging.info("ERROR: could not find pre-computed matrices...re-calculating them which may take awhile")
            A, psi = precompute_terms(W, beta_tilde, name, outdir)


    # get meta values from data
    K = len(mu_vec)
    M = beta_tilde.shape[0]

    # intialize values of chain
    logging.info("initializing p")
    p_init = initialize_p(K)
    logging.info("initializing gamma_k")
    logging.info("initializing C_k")
    C_init, gamma_init = initialize_C_gamma(p_init, mu_vec, sigma_vec, M)

    # calculate sigma_e
    sigma_e = (1-ldsc_h2)/float(N)

    # start profile
    pr = cProfile.Profile()
    pr.enable()

    p_est, C_est, weights = gibbs(p_init, gamma_init, C_init, mu_vec, sigma_vec, sigma_e, W, A, psi, beta_tilde, N, its, f)

    # end Profile
    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()
    print_func(s.getvalue(), f)

    print "Estimate: "
    print p_est

    print "C-Estimate: "
    print C_est

    #print beta_tilde
    BURN = its/4
    weights_est = np.divide(weights, its-BURN)

    weights_true = np.asarray(df['BETA_TRUE'])
    #accuracy = np.corrcoef(weights_true, weights_est)
    accuracy = r2_score(weights_true, weights_est)
    print_func("Accuracy: %.4g" % accuracy, f )

    log_like = log_likelihood(beta_tilde, weights_est, sigma_e, W)
    print_func("Log-like: %.4g" % log_like, f)

    # save results in data-frames
    r_df = {'sigma': sigma_vec, 'p': p_est}
    results_df = pd.DataFrame(data=r_df)
    results_file = os.path.join(outdir, name +'.'+str(seed)+'.results')
    results_df.to_csv(results_file, index=False, sep=' ')

    # save estimated effect sizes
    #df2 = {'weights': weights_est}
    #weights_df = pd.DataFrame(data=df2)
    #weights_file = os.path.join(outdir, name +'.'+str(seed)+'.weights')
    #weights_df.to_csv(weights_file, index=False, sep=' ')

    # add weights to the dataframe
    df['WEIGHTS'] = weights_est
    df.to_csv(gwas_file, index=False, sep=' ')

    f.close()

if __name__== "__main__":
  main()
