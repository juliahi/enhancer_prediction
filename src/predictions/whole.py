
from ..preprocessing import load_vista
from ..shared import *
from ..promoters import promoters


from numpy import transpose, concatenate, array, interp, linspace, hstack


import sys
import os


from ..rf.read_data import read_data, join_and_balance, join_data, cut_rows
#from ..promoters.promoters import *
from ..promoters import read_wiggle
#from ..predictions.predict import find_predictions



def load(histmods, kmers, chrom_name, dist=False):
    target = load_vista.load_enhancers(WIGGLENAME+chrom_name)
    data = read_data(histmods, kmers, target, WIGGLENAME+chrom_name)
    return data


def load_all(histmods, kmers, tissue='', dist=False):
    result = None
    if tissue != '': return load(histmods, kmers, tissue)
    for chrom in range(1,23) + ['X']:
        data = load(histmods, kmers, "chr"+str(chrom)+'.id')

        if result is None: result = data
        else: result = join_data(result, data)
        print chrom, "loaded"
    return result



def load_target(chrom_name, dist=False):
    target_data = load_vista.load_enhancers_with_seq(WIGGLENAME+chrom_name)
    return target_data
    

def load_target_all(tissue=''):
    if tissue != '': return load_target(tissue)
    result = []
    for chrom in range(1,23) + ['X']:
        result += load_target("chr"+str(chrom)+'.id')
    return result


#def load_positions(chrom):
    #return [x.start for x in load_vista.load_enhancers(chrom)]

#def load_ids(chrom):
    #return [(x.start, x.id) for x in load_vista.load_enhancers(chrom)]


