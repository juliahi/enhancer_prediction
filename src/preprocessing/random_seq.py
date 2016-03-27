

from load_vista import load_enhancers, Position, save_enhancers, tissue_types, filter_type
import numpy
import os
import scipy
from scipy import stats
import random
from Bio import SeqIO
from random import shuffle
import fileinput
from ..shared import *

import sys


from fit_distr import *

from overlaps import overlap
from ..rf.load import load_target
import histmods


def gen_randoms(enhancers, enhancers_all, randoms_file, id_start):
    #generate 
    lens = get_lengths(enhancers)
    distr = get_distr(lens)
    
    to_return = []
    e_id = id_start

    record_dict = SeqIO.index(DATAPATH+GENOME_FILE, "fasta")
    l_list = list(distr(size=len(lens)))
    
    print "genome loaded"
    
    for e in enhancers:
        l = l_list.pop()
        while True:
            start = random.randint(0,chr_lens[e.chromosome]-l-1)
            end = start + l - 1
            print l, start, end

            if not overlap(enhancers_all, e.chromosome, start, end):
                seq = record_dict[e.chromosome].seq[start:(end+1)]
                if seq.find('N') == -1:
                    #if histmods.check_position(e.chromosome, start, end):
                        to_return.append(Position(e.chromosome, start, end, e_id, False, [], seq))
                        #print "accepted"
                        break
        e_id+=1
    
    save_enhancers(to_return, randoms_file)
    return to_return


def extend(enhancers, model, outname):
    #extend enhancers sequences to match length of model sequences
    
    enhancers = sorted([x for x in enhancers ], key=lambda e: e.end - e.start)
    lens = get_lengths(model)
    
    m, v, p, r = dist_parameters(lens)

    to_return = []

    l_list = sorted(list(stats.nbinom(r, 1-p).rvs(size=len(enhancers))))
    
    record_dict = SeqIO.index(DATAPATH+GENOME_FILE, "fasta")
    for e in enhancers:
        l = l_list.pop()/2
        start = max(e.start-l,1)
        end = min(e.end+l, chr_lens[e.chromosome]) 
        
        to_return.append(Position(e.chromosome, start, end , e.id+10000, False, [], record_dict[e.chromosome].seq[start:end]))

    save_enhancers(to_return, outname)
    


if  __name__ =='__main__':

    if len(sys.argv) < 4:
        print "USAGE: random_seq.py gen db tissue outname start_id"
        print "OR: random_seq.py extend to_change_db tissue outname model_db tissue"
        sys.exit(1)
    
    outname = sys.argv[4]
    enhancers = load_target(sys.argv[2], sys.argv[3], True)    
    enhancers_all = load_target(sys.argv[2], 'both', True)

    if sys.argv[1] == 'gen':
        gen_randoms(enhancers, enhancers_all, outname, int(sys.argv[5]))
    elif sys.argv[1] == 'extend':
        model = load_target(sys.argv[4], sys.argv[5], True)
        extend(enhancers, model, outname)
    
    
    
    
    
