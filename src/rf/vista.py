
#import matplotlib
#matplotlib.use('Agg')
from ..preprocessing import load_vista
#from boruta import *
from ..shared import *
#import new_train as train
from ..promoters import promoters


from numpy import transpose, concatenate, array, interp, linspace, hstack
#import numpy
#from glob import glob
#from random import sample
#from matplotlib import pyplot as plt
#import pyroc
#import pickle

import sys


from read_data import read_data, join_and_balance


def get_tissue(tissue, dist, not_promoters = False):
    if not dist:
        if tissue == "brain":
            return lambda x:x.is_neural()
        elif tissue == "limb":
            return lambda x:x.is_limb()
        elif tissue == "heart":
            return lambda x:x.is_heart() 
        elif tissue == "positives":
            return lambda x:x.positive
        elif tissue == "vista":
            return lambda x:(not x.is_random())
        elif tissue == "both":
            return lambda x:(x.is_heart() or x.is_neural())
        elif tissue == "all":
            return lambda x: True
    else:
        if tissue == "brain":
            return lambda x:(x.is_neural() and not x.is_heart())
        elif tissue == "heart":
            return lambda x:(x.is_heart() and not x.is_neural())
        elif tissue == "both":
            return lambda x:((x.is_heart() and not x.is_neural()) or (not x.is_heart() and x.is_neural()))        
        elif tissue == "all":
            return lambda x: True
        
    assert False, ("vista option of tissue, dist = ", tissue, dist, "not defined")



def has_tss(enh, tss_list):
    for x in tss_list[enh.chromosome]:
        if enh.start < x[0] < enh.end:
            return True
    return False

def choose_tissue(data, target, tissue, dist, filter_promoters=False):
    fun = get_tissue(tissue, dist)
    tissue_vec = [x for x in range(len(target)) if fun(target[x])]
    if filter_promoters:
        tss = promoters.read_tss()
        prev_len = len(tissue_vec)
        tissue_vec = [x for x in tissue_vec if not has_tss(target[x], tss)]
        print "filtered %d sequences without tss from %d" % (len(tissue_vec), prev_len)
    data = data[[tissue_vec]]
    return data



def load(histmods, kmers, tissue, dist):
    if '_notss' in tissue:
        filter_promoters = True
        tissue = tissue[:-6]
    else:
        filter_promoters = False
    if 'random' in tissue:
            target_data = load_vista.load_enhancers_with_seq(tissue)
            data = read_data(histmods, kmers, target_data, tissue)
            data = choose_tissue(data, target_data, 'all', dist, filter_promoters)
    else:
        target_data = load_vista.load_enhancers_with_seq(VISTA_FILE)
        data = read_data(histmods, kmers, target_data, VISTA_FILE)
        data = choose_tissue(data, target_data, tissue, dist, filter_promoters)
    
    return data

def load_other(database, histmods, kmers, tissue, dist):
    if '_notss' in tissue:
        filter_promoters = True
        tissue = tissue[:-6]
    else:
        filter_promoters = False
    if 'randoms' in tissue:
            target_data = load_vista.load_enhancers_with_seq(tissue+'1500')
            data = read_data(histmods, kmers, target_data, tissue+'1500')
            data = choose_tissue(data, target_data, 'all', dist, filter_promoters)
    else:
        target_data = load_vista.load_enhancers_with_seq(database)
        data = read_data(histmods, kmers, target_data, database)
        data = choose_tissue(data, target_data, tissue, dist, filter_promoters)
    return data




def load_target(tissue, dist):
    if '_notss' in tissue:
        filter_promoters = True
        tissue = tissue[:-6]
    else:
        filter_promoters = False
    if 'random' in tissue:
        #if tissue == 'randoms_brain': tissue = 'randoms_neural'
        target = load_vista.load_enhancers_with_seq(tissue)
        if filter_promoters:
            tss = promoters.read_tss()
            return [x for x in target if not has_tss(x, tss)]
        return target
    else:
        target = load_vista.load_enhancers_with_seq(VISTA_FILE)
        
        fun = get_tissue(tissue, dist)
        
        if filter_promoters:
            tss = promoters.read_tss()
            return [x for x in target if fun(x) and not has_tss(x, tss)]
        return [x for x in target if fun(x)]
    
    
def load_target_other(database, tissue, dist):
    if tissue == 'randoms_both':
        target = load_vista.load_enhancers_with_seq('randoms_heart1500') + load_vista.load_enhancers_with_seq('randoms_brain1500')
    elif 'randoms' in tissue:
        target = load_vista.load_enhancers_with_seq(tissue+'1500')
    else:
        target = load_vista.load_enhancers_with_seq(database)
        if '_notss' in tissue:
            tissue = tissue[:-6]
            tss = promoters.read_tss()
            fun = get_tissue(tissue, dist)
            return [x for x in target if fun(x) and not has_tss(x, tss)]
        fun = get_tissue(tissue, dist)
        target = [x for x in target if fun(x)]
    return target
