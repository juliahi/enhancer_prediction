
from ..preprocessing import load_vista
from ..shared import *
from ..promoters import promoters
from ..promoters import read_wiggle


from numpy import transpose, concatenate, array, interp, linspace, hstack


import sys
import os
from glob import glob

from ..rf.read_data import read_data, join_and_balance, join_data, cut_rows



#def load(histmods, kmers, chrom, dist=False):
    #target_data = load_vista.load_enhancers_with_seq(chrom)
    #data = read_data(histmods, kmers, target_data, chrom)
    #return data


def load_data(histmods, kmers, tissue):
    result = None
    
    target = load_target_dict(tissue)
    
    
    
    for chrom in range(1,23) + ['X']:

        targets_starts = [x.start for x in target['chr'+str(chrom)]]
        
        starts = [x.start for x in load_vista.load_enhancers(WIGGLENAME+"chr"+str(chrom)+'.id')]
        data = read_data(histmods, kmers, starts, WIGGLENAME+'chr'+str(chrom)+'.id')
        #print targets_starts
        #print starts
        
        rows = [i for i in xrange(data.shape[0]) if starts[i] in targets_starts]
        #print rows
        data = cut_rows(data, rows)
        
        if result is None: result = data
        else: result = join_data(result, data)
        print chrom, "loaded"
    return result



def load_target(tissue):
    result = []
    print DATAPATH+tissue+"/*.bed", glob(DATAPATH+tissue+"/*.bed")
    for filename in glob(DATAPATH+tissue+"/*.bed"):
        result += read_wiggle.read_bed(filename)
    return result

def load_target_dict(tissue):
    result = {}
    print DATAPATH+tissue+"/*.bed", glob(DATAPATH+tissue+"/*.bed")
    for filename in glob(DATAPATH+tissue+"/*.bed"):
        result.update(read_wiggle.read_bed_dict(filename))
    
    return result



def load_positions(chrom):
    return [x.start for x in load_vista.load_enhancers(chrom)]

#def load_ids(chrom):
    #return [(x.start, x.id) for x in load_vista.load_enhancers(chrom)]


def predictions_name(histmods, kmers, tissue):
    name = 'predicted'
    if histmods != '-' and histmods != '':
            name += '_'+histmods.split('.')[0]
    if kmers != '-' and kmers != '' and kmers != []:
        if type(kmers)==list:
            name += '_'+kmers[0]
        else:
            name += '_'+str(kmers)
    name += '_' + tissue + '_0.8'
    return name
