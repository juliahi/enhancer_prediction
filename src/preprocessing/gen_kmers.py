from load_vista import load_enhancers, load_enhancers_with_seq, save_enhancers
#from lens import chr_lens
#from bx.bbi.bigwig_file import BigWigFile
from glob import glob
import os
import StringIO
import time
import sys
from Bio import SeqIO

from ..shared import *

def complement(x):
    if x == 'A': return 'T'
    elif x == 'C': return 'G'
    elif x == 'G': return 'C'
    elif x == 'T': return 'A'
    else: 
        print x
        return x

def is_compl(x, y):
    if len(x) != len(y):
        return False
    for i in xrange(len(x)):
        if complement(x[i]) != y[-i-1]: 
            return False
    return True

def get_compl(x): 
    y = ''
    for a in x:
        y = complement(a) + y
    return y


def gen_all(k):
    d = ['']
    alph = ['A', 'C', 'G', 'T']
    for i in xrange(k):
        new_d = []
        for s in d:
            for a in alph:
                new_d.append(s+a)
        d = new_d
        
    #filter reverse complementary
    final = []
    for seq in d:
        to_add=True
        for seq2 in final:
            if is_compl(seq, seq2):
                to_add = False
                break
        if to_add:
            final.append(seq)
    #print final
    print "generated %d %d-mers" % (len(final), k)
    return final


def change(kmers, s, val):
    #change value in kmers dictionary for sequence s (or complementary sequence) by val
    #return new value
    if kmers.has_key(s):
        kmers[s] += val
        return kmers[s]
    kompl = get_compl(s)
    if kmers.has_key(kompl):
        kmers[kompl] += val
        return kmers[kompl]
    print "error: %s not in kmers dictionary" % s
    
    
def count(keys, k, seq):
    kmers = dict.fromkeys(keys, 0)
    for i in xrange(0, len(seq) - k+1):   #zmiana 1.12.15 - dodanie + 1
            s = str(seq[i:(i+k)].upper())
            change(kmers, s, 1)
    if sum(kmers.values()) != len(seq)-k+1:
        print "missing kmers", k, sum(kmers.values()), len(seq)
        
    return kmers



def count_kmers(enhancers_file, k):
    print "Counting %d-mers for file %s" % (k, enhancers_file)
    keys = gen_all(k) 
    enhancers=load_enhancers_with_seq(enhancers_file)
    
    output = open(DATAPATH+enhancers_file+'.%dmers' % k, "w")
    for x in keys:
        output.write("\t%s" %x)
    output.write('\n')
    for enh in enhancers:
        output.write("%s\t" % enh.id)
        assert len(enh.seq) > 0,             "no sequence for %d"%enh.id

        kmers = count(keys, k, enh.seq)
        for key in keys:
                output.write("%f\t" % (float(kmers[key])/len(enh.seq)))
        output.write('\n')
    output.close()



def count_kmers_maxes(enhancers_file, k, nmaxes, window=200):
    print "%d-MERY z %d maksami w oknie %d w pliku %s" % (k, nmaxes, window, enhancers_file)
    keys = gen_all(k) 
    enhancers=load_enhancers_with_seq(enhancers_file)
    
    output = open(DATAPATH+enhancers_file+'.%dmaks%d_%d' % (k, nmaxes, window), "w")
    #header
    
    for x in keys:
        for maks in xrange(nmaxes):
            output.write("\t%s-%dmaks" % (x, maks))
    output.write('\n')
    #count for every sequence
    for enh in enhancers:
        output.write("%s\t" % enh.id)
        assert len(enh.seq) > 0,             "no sequence for %d"%enh.id
        
        n = len(enh.seq)
        kmers = count(keys, k, enh.seq[:window])
        maxes = dict([(key, [v]) for key,v in kmers.items()])
        for i in xrange(0, n - window - k, 1):
            sold = str(enh.seq[i:(i+k)].upper())
            snew = str(enh.seq[(i+window):(i+window+k)].upper())
            
            x = change(kmers, sold, -1)
            y = change(kmers, snew, 1)
            
            if sold in maxes:
                maxes[sold].append(x)
            elif get_compl(sold) in maxes:
                maxes[get_compl(sold)].append(x)
            
            if snew in maxes:
                maxes[snew].append(y)
            elif get_compl(snew) in maxes:
                maxes[get_compl(snew)].append(y)
            
            
        for key in maxes:
                sortedmax = sorted(maxes[key] + [0]*nmaxes, reverse=True)[:nmaxes]
                for maks in xrange(nmaxes):
                    output.write("%f\t" % (float(sortedmax[maks])/window))
        output.write('\n')
    output.close()


    
    
if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print "USAGE: gen_kmers.py input_filename K [n_maxes window] "
        print "suggested: window=100 for FANTOM, window = 300 for VISTA"
        sys.exit(1)
    f = sys.argv[1]
    k = int(sys.argv[2])
    if len(sys.argv) == 3:
        count_kmers(f, k)
    else:
        count_kmers_maxes(f, k, int(sys.argv[3]), int(sys.argv[4]))

