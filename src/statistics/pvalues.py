



from ..shared import *

import roc_compare
import pandas

import itertools
import numpy
from rpy2 import robjects
from ..rf import read_data, load
from ..rf import new_train as train
import math
from scipy.stats.stats import pearsonr

import summary



def compare_pvalues(db, tissue):
    predictions = {}
    options = [('h1hesc.list', '-'), 
                            ('h1hesc.list', '4mers'), 
                            ('-', '4mers'), 
                            ('tier12.list', '-'), 
                            ('tier12.list', '4mers'), 
                            ]
    for histmods, kmers in options:
        
        datapos = load.load_data(db, histmods, kmers, tissue, dist=True )
        dataneg = load.load_data(db, histmods, kmers, 'randoms_'+tissue, dist=True )
        
        print datapos.shape, dataneg.shape
        npos = datapos.shape[0]
        total = 2*npos
        
        target = [1]*npos + [0]*npos
        data = read_data.join_data(datapos, dataneg)
        data, names = read_data.my_transpose(data)
        scores = train.predict_cv(data, N_trees, target)
        
        predictions[(histmods, kmers)] = scores
    
    
    with open(RESULTSPATH+"pvalues_cROC.csv", "a") as f:
        for opt1, opt2 in itertools.combinations(options,2):
            if set('-').issuperset(set(opt1) & set(opt2)):
                continue
            pred1 = predictions[opt1]
            pred2 = predictions[opt2]
            
            auc1 = train.auc2(pred1[0], pred1[1])
            auc2 = train.auc2(pred2[0], pred2[1])
            
            if (auc1 > auc2):
                auc2, auc1 = auc1, auc2
                pred2, pred1 = pred1, pred2
                opt2, opt1 = opt1, opt2
                
            #r = robjects.r.cor(pred1[0] + pred1[1],pred2[0] + pred2[1])[0]
            r = pearsonr(pred1[0] + pred1[1],pred2[0] + pred2[1])[0]
            
            pvalue = roc_compare.cROC(auc1,auc2,npos,total,r)
            #print opt1, opt2
            #if pvalue < 0.05:
                #print '-------------------', auc1, auc2, r, pvalue
            #else:
                #print auc1, auc2, r, pvalue
            
            
            
            pvalue2 = roc_compare.diffAUC(auc2,auc1,npos,total)
            conf = ''
            if pvalue2 < 0.05:
                conf += '*'
            if pvalue2 < 0.005:
                conf += '*'
            if pvalue2 < 0.0005:
                conf += '*'
            #conf = int(math.log10(pvalue2/5)-1)
            f.write("\t".join([db, tissue, opt1[0], opt1[1], opt2[0], opt2[1]] + map(str, [auc1, auc2, r, pvalue, pvalue2, conf])) + "\n")

def compare_pvalues_mean(db, tissue):
    options = [('h1hesc.list', '-'), 
                            ('h1hesc.list', "['4mers']"), 
                            ('-', "['4mers']"), 
                            ('tier12.list', '-'), 
                            ('tier12.list', "['4mers']"), 
                            ('ALL', "-"), ('ALL', "['4mers']"), ('No', "['4mers']"), 
                            ]
    mean_aucs = summary.read_summary(db=db, tissue=tissue)
    print mean_aucs
    
    with open(RESULTSPATH+"pvalues_diffAUC.csv", "a") as f:
        for opt1, opt2 in itertools.combinations(options,2):
            if set('-').issuperset(set(opt1) & set(opt2)):
                continue
            
            
            try:
                auc1 = mean_aucs.loc[opt1[1],opt1[0]]['AUC']['mean']
                auc2 = mean_aucs.loc[opt2[1],opt2[0]]['AUC']['mean']
            except:
                continue
            print auc1, auc2
            if (auc1 > auc2):
                auc2, auc1 = auc1, auc2
                opt2, opt1 = opt1, opt2
                
            npos = mean_aucs.loc[opt1[1],opt1[0]]['POS_size']
            nneg = mean_aucs.loc[opt1[1],opt1[0]]['NEG_size']
            pvalue = roc_compare.diffAUC(auc2,auc1,npos,npos+nneg)
            conf = ''
            if pvalue < 0.05:
                conf += '*'
            if pvalue < 0.005:
                conf += '*'
            if pvalue < 0.0005:
                conf += '*'
            #conf = int(math.log10(pvalue2/5)-1)
            f.write("\t".join([db, tissue, opt1[0], opt1[1], opt2[0], opt2[1]] + map(str, [auc1, auc2, pvalue, conf])) + "\n")
        f.write('\n')


if __name__ == '__main__':
    for db in ['vista', 'fantom']:
        compare_pvalues_mean(db, 'heart')
        compare_pvalues_mean(db, 'brain')
        compare_pvalues_mean(db, 'both')