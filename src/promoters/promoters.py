
from ..shared import *
import read_wiggle
from ..rf.read_data import read_data, join_and_balance, join_data, cut_rows
from ..preprocessing.load_vista import load_enhancers

import pickle
import numpy
import os
import sys

def read_tss():
    promoters = {}
    with open(TSS_FILE, 'r') as f:
        next(f)
        for line in f:
            prom = line.strip().split('\t')
            if len(prom) < 6:
                continue
            if prom[5] == '+' or prom[5] == "1":
                pos = (int(prom[1]), 1)
            else:
                pos = (int(prom[2]), -1)
            
            if prom[0] == 'MT': prom[0] = 'M'
            if prom[0][:3] != 'chr':
                prom[0] = 'chr'+prom[0]
            if prom[0] in promoters:
                promoters[prom[0]].append(pos)
            else:
                promoters[prom[0]] = [pos]
    #print promoters.keys()
    return promoters




def promoters_name(histmods, kmers, tissue):
    name = 'predicted_promoters'
    if histmods != '-' and histmods != '':
            name += '_'+histmods.split('.')[0]
    if kmers != '-' and kmers != '':
            name += '_'+kmers
    name += '_' + tissue + '_0.8'
    return name




def load_target(tissue):
    result = []
    for chrom in range(1,23) + ['X']:
        if os.path.exists(DATAPATH+tissue+"/chr"+str(chrom)+"_sorted.bed"):
            result += read_wiggle.read_bed(DATAPATH+tissue+"/chr"+str(chrom)+'_sorted.bed')
    return result
        

#def load_all(length = 1500):
    ##returns list of regions of length length ending at tss (all promoters)
    #with open(TSS_FILE, 'r') as f:
        #next(f)
        #promoters = []
        #for line in f:
            #prom = line.strip().split('\t')
            #if len(prom) < 6:
                #continue
            
            #if prom[0] == 'MT': prom[0] = 'M'
            #if len(prom[0]) <= 3:
                #prom[0] = 'chr'+prom[0]
            #if prom[5] == '1' or prom[5] == '+':
                #pos = Region(prom[0], int(prom[1])-length, int(prom[1]), prom[5] )
            #else:
                #pos = Region(prom[0], int(prom[2]), int(prom[2])+length, prom[5] )
            #promoters.append(pos)
    #return promoters



#def load_dict(tissue):
    #if tissue == 'all':
        #return load_all()
    #else: #load predictions from filename = tissue
        #with open(DATAPATH+tissue+'_sorted.bed', 'r') as f:
            #promoters = {}
            #for line in f:
                #prom = line.strip().split('\t')
                #if len(prom) < 6:
                    #continue
                #chrom = prom[0]
                #pos = Region(prom[0], int(prom[1]), int(prom[2])-1, prom[5] )
                #if chrom in promoters:
                    #promoters[chrom].append(pos)
                #else:
                    #promoters[chrom] = [pos]
        #return promoters


#def load_michal(histmods, kmers, tissue):
    #targets = load_dict(tissue)
    
    #try:
        #return pickle.load(open(DATAPATH+'load_promoters_%s_%s_%s'%(tissue, histmods, kmers), 'r'))
    #except:
    
        #data = None
        #for chrom, promoters in targets.iteritems():
            #data2 = michal_read_data.read_data(chrom, histmods, kmers, [x.start for x in promoters])
            #print chrom, data2.shape
            #if data != None:
                #data = numpy.append(data, data2)
            #else:
                #data = data2
        
        #pickle.dump(data, open(DATAPATH+'load_promoters_'+tissue, 'w'))
        #return data
        
        
#def load(histmods, kmers, tissue):    
    #result = None
    
    
    #for chrom in range(1,23) + ['X']:
        ##print DATAPATH+tissue+"/chr"+str(chrom)+'_sorted.bed', os.path.exists(DATAPATH+tissue+"/chr"+str(chrom)+'_sorted.bed')
        #if not os.path.exists(DATAPATH+tissue+"/chr"+str(chrom)+'_sorted.bed'): continue
        #print "loading", chrom
        #targets = read_wiggle.read_bed(DATAPATH+tissue+"/chr"+str(chrom)+'_sorted.bed')
        #targets_starts = [x.start for x in targets]
        
        #starts = [x.start for x in load_enhancers(WIGGLENAME+"chr"+str(chrom)+'.id')]
        #data = read_data(histmods, kmers, starts, WIGGLENAME+'chr'+str(chrom)+'.id')
        ##print targets_starts
        ##print starts
        
        #rows = [i for i in xrange(data.shape[0]) if starts[i] in targets_starts]
        #print rows
        #data = cut_rows(data, rows)
        
        #if result is None: result = data
        #else: result = join_data(result, data)
        #print chrom, "loaded"
    #return result
    
    
    
