

from random import sample
import numpy

from ..shared import *
import vista
import fantom
from ..promoters import promoters
from ..predictions import predictions, whole
import os
import pickle
import read_data
import train_two_step 

from ..promoters import read_wiggle
from ..preprocessing.load_vista import load_enhancers
    

def load_target(database, tissue, dist = True):
    if 'balanced' in tissue:
        data1 = load_target(database, tissue.replace('balanced', 'heart'), dist)
        data2 = load_target(database, tissue.replace('balanced', 'brain'), dist)
        n = min(len(data1), len(data2))
        idx1 = numpy.random.choice( len(data1), n, replace=False)
        idx2 = numpy.random.choice( len(data2), n, replace=False)
        return [data1[x] for x in idx1] +  [data2[x] for x in idx2] 
    if 'both' in tissue and database != 'promoters':
        data1 = load_target(database, tissue.replace('both', 'heart'), dist)
        data2 = load_target(database, tissue.replace('both', 'brain'), dist)
        return data1 + data2
    
    if database == 'vista':
        data = vista.load_target(tissue, dist)
    if database == "vista1500":
        data = vista.load_target_other(database, tissue, dist)
    if database == 'fantom':
        data = fantom.load_target(tissue, True, dist)
    if database == 'fantom_long':
        data = fantom.load_target(tissue, False, dist)
    if database == 'promoters':
        data = promoters.load_target(tissue)
    if database == 'predictions':
        data = predictions.load_target(tissue)
    if database == 'whole':
        data = whole.load_target_all(tissue)
        print tissue
        
    print database, tissue
    return data
    


n = 50 #for balanced sets

def load_data(database, histmods, kmers, tissue, dist = True):
    print tissue, histmods, kmers

    if 'balanced' in tissue:
        data1 = load_data(database, histmods, kmers, tissue.replace('balanced', 'heart'), dist)
        data2 = load_data(database, histmods, kmers, tissue.replace('balanced', 'brain'), dist)
        n = min(data1.shape[0], data2.shape[0])
        return read_data.join_data(train_two_step.choose(data1, n), train_two_step.choose(data2, n))
    if 'both' in tissue and database != 'promoters':
        data1 = load_data(database, histmods, kmers, tissue.replace('both', 'heart'), dist)
        data2 = load_data(database, histmods, kmers, tissue.replace('both', 'brain'), dist)
        return read_data.join_data(data1, data2)
    
    if kmers == '-' or kmers == '': kmers = []
    if type(kmers) != list:
        kmers = [kmers]
    if histmods == '-': histmods = ''
    
    path = database+'_'+histmods+'_'+str(kmers)+'_'+tissue+'.pickle'
    path = RESULTSPATH+'datasets/'+path.replace('/', '_')
    if os.path.exists(path):
        data = pickle.load(open(path))
        print 'loading saved', path
        return data
    
    
    if database == 'vista':
        data = vista.load(histmods, kmers, tissue, dist)
    elif database == "vista1500":
        data = vista.load_other(database, histmods, kmers, tissue, dist)
    elif database == 'fantom':
        data = fantom.load(histmods, kmers, tissue, True, dist)
    elif database == 'fantom_long':
        data = fantom.load(histmods, kmers, tissue, False, dist)
    
    elif database == 'promoters':
        result = None
        for chrom in range(1,23) + ['X']:
        #print DATAPATH+tissue+"/chr"+str(chrom)+'_sorted.bed', os.path.exists(DATAPATH+tissue+"/chr"+str(chrom)+'_sorted.bed')
            if not os.path.exists(DATAPATH+tissue+"/chr"+str(chrom)+'_sorted.bed'): continue
            print "loading", chrom
            targets = read_wiggle.read_bed(DATAPATH+tissue+"/chr"+str(chrom)+'_sorted.bed')
            targets_starts = [x.start for x in targets]
        
            starts = [x.start for x in load_enhancers(WIGGLENAME+"chr"+str(chrom)+'.id')]
            if 'heart' in tissue: t = 'heart'
            else: t = 'brain'
            data = load_data('whole', histmods, kmers, "chr"+str(chrom)+'.id')
            #print targets_starts
            #print starts
            print len(targets_starts)
        
            rows = [i for i in xrange(data.shape[0]) if starts[i] in targets_starts]
            #print rows
        
            data = read_data.cut_rows(data, rows)
        
            if result is None: result = data
            else: result = read_data.join_data(result, data)
            print chrom, "loaded"
        data=result

    elif database == 'predictions':
        print tissue, histmods, kmers
        data = predictions.load_data(histmods, kmers, tissue)
        print data    
    elif database == 'whole':
        print tissue, histmods, kmers
        data = whole.load_all(histmods, kmers, tissue)
    else:
        print database, tissue, "not defined"
 
    print data.shape
    pickle.dump(data, open(path, 'w+'))
    return data







    
