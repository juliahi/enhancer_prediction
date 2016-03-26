

from ..shared import *

from random import sample
from glob import glob
import numpy


def read_data(filename, kmers_filenames, sequences, name):
    
    names = []
    knames = []
    if kmers_filenames == '':
        useKmers = False
    else:
        useKmers = True
        kmers = []
        for kmers_filename in kmers_filenames:
            kmers += [open(DATAPATH+name+"."+kmers_filename)]
            knames += kmers[-1].readline().rstrip().split('\t')[1:]
    if filename != '':
        directories = []
        with open(DATAPATH+filename) as f: directories += map(lambda x:x.strip(), f.readlines())
        filenames = []
        for directory in directories:
            filenames += sorted(glob(DATAPATH+directory + "/*.means." + name.split('/')[-1]))
        for filename in filenames:
            names.append(get_name(filename))
        #print filenames
    else:
        directories = []
        filenames = []
    
    
    if USEGC:
        names += knames + ["GC-content"]
    else:
        names += knames
    data = numpy.zeros(shape=(len(sequences),), dtype={'names': names, 'formats': ['float']*len(names)})
    #load histone mods
    for filename in filenames:
            with open(filename) as f:
                row = []
                for line in f.readlines():
                    assert len(line.strip().split()) > 1, "file {0} corrupted! Line not finished: {1}".format(filename, line)
                    row.append(line.strip().split()[1])
                
                assert len(sequences) == len(row), "file {0} corrupted! Length {1} not {2}".format(filename, len(row), len(sequences))
                data[get_name(filename)] = row      
    
    #load kmers
    print filenames, knames
    j = len(filenames)
    if useKmers:
        for i in xrange(len(sequences)):
            x = []
            
            for fileid in xrange(len(kmers)):
                line = kmers[fileid].readline()
                x+= map(float, line.split()[1:])
            #row = [a[0] + a[1] for a in zip(data[i], [0]*j+x+[0])]
            row = tuple(data[i])[:j]+ tuple(x) #zmiana!
            data[i] = tuple(row)
            
    
    #add lengths and GC-content
    #lengths = []
    if USEGC:
        gc = []
        for e in sequences:
            #    lengths.append(len(e.seq))
            gc.append(float(e.seq.count("G") + e.seq.count("C") + e.seq.count("g") + e.seq.count("c"))/len(e.seq))
        #data["lengths"] = lengths
        data["GC-content"] = gc
    return data


def cut_rows(data, rows):
    return data[[rows]]

def join_data(data1, data2):
    return numpy.append(data1, data2)

def join_and_balance(data1, data2, balance=True):
        # zrownowazenie wielkosci grup: 
    #print "join and balance", data1.shape, data2.shape
    len_pos = data1.shape[0]
    len_neg = data2.shape[0]

    if balance and len_pos != len_neg:
        len_neg2 = min(len_neg, len_pos ) 
        len_pos2 = min(len_neg, len_pos )  
        pos_i = sample(range(len_pos), len_pos2)
        neg_i = sample(range(len_neg), len_neg2)
        data = cut_rows(data1, pos_i)
        data = numpy.append(data, cut_rows(data2, neg_i))
    else:
        len_neg2 = len_neg
        len_pos2 = len_pos
        data = join_data(data1, data2)
    
    
    print "data sizes %d->%d, %d->%d" % (len_pos, len_pos2, len_neg, len_neg2)
    target = [1]*len_pos2 + [0]*len_neg2
    return data, target


def my_transpose(data):
    #print data.shape[0], len(data.dtype.names), data[0]
    return data.view(float).reshape(data.shape[0], len(data.dtype.names)), data.dtype.names
    
    
    #array = numpy.empty((data.shape[0], len(data.dtype.names)))
    #for x in xrange(data.shape[0]):
        #for i in xrange(len(array[x])):
            #array[x][i] = data[x][i]
    #return array, data.dtype.names
    
