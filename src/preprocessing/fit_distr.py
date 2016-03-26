import scipy
from scipy import stats
from ..shared import *


def get_distr_fun(lens):
    _, _, p, r = dist_parameters(lens)
    return stats.nbinom(r, 1-p).pmf

def get_distr(lens):
    _, _, p, r = dist_parameters(lens)
    return stats.nbinom(r, 1-p).rvs


def dist_parameters(lens):
    #compute Negative binomial distribution parameters

    m = scipy.mean(lens)
    v = scipy.var(lens)
    p = (v-m)/v
    r = m*(1-p)/p
    
    return m, v, p, r


def get_lengths(seqs):
    return  map(lambda x:x.end - x.start + 1, seqs)

def print_stats(seqs):
    lens = get_lengths(seqs)
    m, v, p, r = dist_parameters(lens)
    print "mean\tvariance\tmedian\tr\tp"
    print "\t".join(map(str, [m, v, scipy.median(lens), r, p]))
    
    return lens, stats.nbinom(r, 1-p).pmf
