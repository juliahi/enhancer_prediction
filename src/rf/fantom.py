# -*- coding: utf-8 -*-

#import matplotlib
#matplotlib.use('Agg')


#from boruta import *
from ..shared import *
from ..preprocessing import load_vista 


import os.path
#import new_train as train


from read_data import read_data, join_and_balance


def remove_samples(data, target_data, files):
    #remove samples from data with target_data ids listed in files
    not_remove = []
    lista = []
    for f in files:
        if os.path.exists(DATAPATH+f):
            lista += [int(x) for x in open(DATAPATH+f, "r").readlines() ]
    i = 0
    for e in target_data:
        if not (e.id ) in lista:
            not_remove.append(i)
        i += 1
    return data[[not_remove]]


def get_filename(tissue, short):
    if "random" in tissue:
        filename = "FANTOM_" + tissue.split('_')[1] + "_randoms"
        #if short:
        #    filename += ".short"
        if not short:
            filename += ".long"
    else:
        filename = "FANTOM_" + tissue
        if not short:
            filename += ".long"
    return filename



def load_target(tissue, short, dist):
    if tissue == 'both':
        tissues = ['heart', 'brain']
    if tissue == 'randoms_both':
        tissues = ['randoms_heart', 'randoms_brain']
    else:
        tissues = [tissue]
    
    result = []
    for tissue in tissues:
        filename = get_filename(tissue, short)
        
        target_data = load_vista.load_enhancers_with_seq(filename)
        #24.07.2015 - zmiana --> endpoint w plikach FANTOM nie jest zawarty w sekwencji (inaczej niz vista)
        for e in target_data:
            e.end = e.end-1
        
        files_remove = [filename+".remove"]
        if dist:
            if tissue == "heart":
                files_remove.append( "FANTOM_heart_brain.remove" )
            elif tissue == "brain":
                files_remove.append("FANTOM_brain_heart.remove")
        
        not_remove = []
        lista = []
        if 'randoms' not in tissue:   #no .remove for randoms files YET
            for f in files_remove:
                lista += [int(x) for x in open(DATAPATH+f, "r").readlines() ]
        
            for e in target_data:
                if not (e.id ) in lista:
                    result.append(e)
        else:
            result = target_data
    return result


def load(histmods, kmers, tissue, short, dist):
    if 'both' in tissue:
        if 'randoms' in tissue:
            randoms = 'randoms_'
        else: randoms = ''
        filename = get_filename(randoms +'heart', short)
        filename2 = get_filename(randoms +'brain', short)
        target_data = load_vista.load_enhancers_with_seq(filename)
        target_data2 = load_vista.load_enhancers_with_seq(filename2)
        files_remove = [filename+".remove"]
        files_remove2 = [filename2+".remove"]
        if dist:
            files_remove.append( "FANTOM_heart_brain.remove" )
            files_remove2.append("FANTOM_brain_heart.remove")
        
        data = read_data(histmods, kmers, target_data, filename)
        data = remove_samples(data, target_data, files_remove)
        data2 = read_data(histmods, kmers, target_data2, filename2)
        data2 = remove_samples(data2, target_data2, files_remove2)
        return join_and_balance(data, data2, False)[0]
    
    else:
        filename = get_filename(tissue, short)
        target_data = load_vista.load_enhancers_with_seq(filename)
        files_remove = [filename+".remove"]
        if dist:
            if tissue == "heart":
                files_remove.append( "FANTOM_heart_brain.remove" )
            elif tissue == "brain":
                files_remove.append("FANTOM_brain_heart.remove")
        print target_data
        data = read_data(histmods, kmers, target_data, filename)
        data = remove_samples(data, target_data, files_remove)
        return data


