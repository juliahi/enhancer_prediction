

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
                    if histmods.check_position(e.chromosome, start, end):
                        to_return.append(Position(e.chromosome, start, end, e_id, False, [], seq))
                        print "accepted"
                        break
        e_id+=1
    
    save_enhancers(to_return, randoms_file)
    return to_return


#def gen_randoms(enhancers_file):
    #enhancers = load_enhancers(enhancers_file)
    #ids = 100001
    ##for t in tissue_types:
    #for t in ["positives"]:
        #enh_filt = filter_type(enhancers, t)
        #lens = map(lambda x:x.end - x.start, enh_filt)
        #randoms, fun = gen_random(enh_filt, lens, DATAPATH + ("randoms_%s"%t), ids)
        #ids += 10000
        ##lens2 = map(lambda x:x.end - x.start, randoms)
        ##plot_lengths(lens, lens2, fun, DATAPATH + ("img/lengths_%s.pdf"%t))
    
#def gen_fantom2(record_dict, enhancers, lens, randoms_file, id_start):
    #lens2 = map(lambda x:x.end - x.start, enhancers)
    
    #m, v, p, r = dist_parameters(lens)

    #m2, v2, p2, r2 = dist_parameters(lens2)
    
    #to_return = []
    #to_return2 = []
    #e_id = id_start

    #print "generating lengths"
    ##generujemy dwa razy tyle ile trzeba
    #l_list = sorted(list(stats.nbinom(r, 1-p).rvs(size=2*len(lens2))))
    #l_list2 = sorted(list(stats.nbinom(r2, 1-p2).rvs(size=2*len(lens2))))
    #print "generated lengths"
    
    #shuffle(enhancers) 
    #for e in enhancers + enhancers: 
        #l = l_list.pop()
        #l2 = l_list2.pop()
        #while True:
            #start = random.randint(0,chr_lens[e.chromosome]-l-1)
            #seq = record_dict[e.chromosome].seq[start:(start+l)]
            #if seq.find('N') == -1:
                #to_return.append(Position(e.chromosome, start, start + l, e_id, False, [], seq))
                #x = (l-l2)/2
                #to_return2.append(Position(e.chromosome, start+x, start + x+l2, e_id, False, [], seq[x:(x+l2)]))
                #break
        #e_id+=1
    #print "enhancers generated"
    #save_enhancers(to_return, randoms_file)
    #save_enhancers(to_return2, randoms_file+".short")
    #return to_return, stats.nbinom(r, 1-p).pmf


#def gen_randoms_fantom():
    #ids = 100001
    #vista = load_enhancers("new_vista_human")
    #vista_lens = map(lambda x:x.end - x.start, filtered)
    #record_dict = SeqIO.index(DATAPATH+"female.hg19.fa", "fasta")
    #for enhancers_file, filtr in [("heart", lambda x: x.is_heart()), ("brain", lambda x: x.is_neural()), ("permissive", lambda x: True)]:
        #filtered = [x for x in vista if filtr(x)]
        #enhancers = load_enhancers("FANTOM_"+enhancers_file, enhancers_file)
        ##enh_filt = filter_type(enhancers, t)
        ##lens = map(lambda x:x.end - x.start, enhancers)
        #randoms, fun = gen_fantom2(record_dict, enhancers, vista_lens, DATAPATH + ("FANTOM_%s_randoms"%enhancers_file), ids)

        ##diffs = [(l-l2)/2 for l, l2 in zip(vista_lens, lens)]
        #ids += 10000
        ##lens2 = map(lambda x:x.end - x.start, randoms)
        #plot_lengths(vista_lens, randoms, filename=DATAPATH + ("img/lengths_FANTOM_%s.pdf"%enhancers_file))

  

def extend(enhancers, model, outname):
    #extend enhancers sequences to match length of model sequences
    
    enhancers = sorted([x for x in enhancers ], key=lambda e: e.end - e.start)
    lens = get_lengths(model)
    
    m, v, p, r = dist_parameters(lens)

    
    to_return = []

    l_list = sorted(list(stats.nbinom(r, 1-p).rvs(size=len(enhancers))))
    #l_list2 = sorted(list(stats.nbinom(r2, 1-p2).rvs(size=len(lens2))))
    
    
    record_dict = SeqIO.index(DATAPATH+GENOME_FILE, "fasta")
    for e in enhancers:
        l = l_list.pop()/2
        #while True:
        #if record_dict[e.chromosome].seq[start:(start+l)].find('N') == -1:
        start = max(e.start-l,1)
        end = min(e.end+l, chr_lens[e.chromosome]) 
        
        to_return.append(Position(e.chromosome, start, end , e.id+10000, False, [], record_dict[e.chromosome].seq[start:end]))

    save_enhancers(to_return, outname)
    
    
    
#def add_id(name, start):
    ##dodaj kolejne id do pliku fasta z FANTOMa
    #for line in fileinput.input(name, inplace=True):
        #if line[0] == ">":
            #print "%s:%d" % (line[:-1], start)
            #start += 1
        #elif not (line[0] == ":" or line[0] == "/n"):
            #print line[:-1]
    
##def add_id2(name, start):
    ###dodaj kolejne id do pliku fasta z FANTOMa
    ##for line in fileinput.input(name, inplace=True):
        ##if line.startswith("chr"):
            ##print ">%s" % (line[:-1])
        ##else:
            ##print line[:-1]
    



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
    
    
    
    
    
