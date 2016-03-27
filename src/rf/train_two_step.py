# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg', warn=False)


from boruta import *
from ..shared import *
from ..preprocessing import load_vista


from sklearn import svm, grid_search
from sklearn import cross_validation
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree 
from numpy import transpose, concatenate, array, interp, linspace, hstack
import numpy
from glob import glob
from random import sample
from matplotlib import pyplot as plt

import pyroc
import pickle
import sys
import os


from read_data import read_data, join_and_balance, my_transpose, cut_rows
from new_train import cut_at, train_random_forest, get_random_forest, auc

def join(datap1, datan1a, datan1b, datap2, datan2):
    #returns: datap1+datan1a, datap2+datan1b, datap2+datan2 with n=min(all lengths) samples
    n = min(datap1.shape[0], datap2.shape[0], datan1a.shape[0], datan1b.shape[0], datan2.shape[0])
    
    pos_i = sample(range(min(datap1.shape[0], datap2.shape[0])), n)
    neg_i = sample(range(datan1a.shape[0]), n)
    neg_i2 = sample(range(datan2.shape[0]), n)
    
    #sys.stderr.write('%d %d %d %d' %(datap1.shape[0], datap2.shape[0], datan1.shape[0], datan2.shape[0]))
    
    data1 = cut_rows(datap1, pos_i)
    data1 = numpy.append(data1, cut_rows(datan1a, neg_i))
    
    data2p = cut_rows(datap2, pos_i)
    data2 = numpy.append(data2p, cut_rows(datan2, neg_i2))
    
    data1b = numpy.append(data2p, cut_rows(datan1b, neg_i))
    print "data sizes", (data1.shape, data1b.shape, data2.shape )
    #target = [1]*len_pos2 + [0]*len_neg2
    
    
    return my_transpose(data1)[0], my_transpose(data1b)[0], my_transpose(data2)[0]

def combine_predictions(pred, pred2):
    def fun(a,b):
        #if a < 0.5 or b < 0.5:
        #    return 0
        return a*b       
    return [(1-fun(pred[i][1], pred2[i][1]), fun(pred[i][1], pred2[i][1])) for i in xrange(len(pred))]





def predict_cv(pos1, neg1a, neg1b, pos2, neg2):
    #pos1 and pos2 are ordered, each sequence pos1[i] corresponds to i-th sequence in pos2[i]
    print pos1.shape[0], neg1a.shape[0], neg1b.shape[0], pos2.shape[0], neg2.shape[0]

    len_n = min(pos1.shape[0], pos2.shape[0], neg1a.shape[0],  neg1b.shape[0], neg2.shape[0])
    target = [1]*len_n + [0]*len_n
    
    data1a, data1b, data2 = join(pos1, neg1a, neg1b, pos2, neg2)
    
    kf = cross_validation.StratifiedKFold(target, cv_folds)
    predictions_pos = [None]*len_n
    predictions_neg = [None]*len_n

    for train_index, test_index in kf:
        data_train = data1a[[train_index]]
        data_test = data1a[[test_index]]
        data_train2 = data2[[train_index]]
        data_test2 = data1b[[test_index]]
        
        target_train = [target[index] for index in train_index]
        target_test = [target[index] for index in test_index]
        predicted = train_random_forest(data_train, target_train, data_test, N_trees)
        predicted2 = train_random_forest(data_train2, target_train, data_test2, N_trees)
        predicted3 = combine_predictions(predicted, predicted2)
        #predicted3 = predicted   #one step!
        
        for x in xrange(len(test_index)):
            if target[test_index[x]] == 1:
                predictions_pos[test_index[x]] = predicted3[x][1]
                #print '\t'.join(map(str, ['pos', x, test_index[x], target[test_index[x]], predicted[x,1], predicted2[x,1], predicted3[x][1]]))
            else:
                #print '\t'.join(map(str, ['neg', x, test_index[x], target[test_index[x]], predicted[x,1], predicted2[x,1], predicted3[x][1]]))
                predictions_neg[test_index[x]-len_n] = predicted3[x][1]

    return predictions_pos, predictions_neg



def predict_cv_bothpromoters(pos1, neg1a, neg1b, pos2a, neg2a, pos2b, neg2b):
    #pos1 and pos2 are ordered, each sequence pos1[i] corresponds to i-th sequence in pos2[i]
    print pos1.shape[0], pos2a.shape[0], neg1a.shape[0], neg2a.shape[0]
    len_n = min(pos1.shape[0], pos2a.shape[0], neg1a.shape[0], neg2a.shape[0])
    target = [1]*len_n + [0]*len_n
    
    data1a, data1b, data2 = join(pos1, neg1a, neg1b, pos2a, neg2a)
    
    kf = cross_validation.StratifiedKFold(target, cv_folds)
    predictions_pos = [None]*len_n
    predictions_neg = [None]*len_n
    
    

    for train_index, test_index in kf:
        data_train = data1a[[train_index]]
        data_test = data1a[[test_index]]
        data_test2 = data1b[[test_index]]
        
        target_train = [target[index] for index in train_index]
        target_test = [target[index] for index in test_index]
        
        npos = sum([1 for i in train_index if target[i] == 1])
        nneg = sum([1 for i in train_index if target[i] == 0])
        print npos, nneg, pos2b.shape[0], neg2b.shape[0]
        
        
        if npos > pos2b.shape[0] or nneg > neg2b.shape[0]:
            nsamples = min(pos2b.shape[0], neg2b.shape[0])
            new_train_index =  sample(train_index, nsamples)
            data_train2 = data2[[new_train_index]]
            data_train2b, _ = join_and_balance(pos2b[[sample(range(pos2b.shape[0]), nsamples)]], neg2b[[sample(range(neg2b.shape[0]), nsamples)]], balance=False)
            target_train2 =  [target[index] for index in new_train_index] +[1]*nsamples+[0]*nsamples
        else:
            data_train2 = data2[[train_index]]
            data_train2b, _ = join_and_balance(pos2b[[sample(range(pos2b.shape[0]), npos)]], neg2b[[sample(range(neg2b.shape[0]), nneg)]], balance=False)
            target_train2 = target_train +[1]*npos+[0]*nneg
        
        data_train2b = my_transpose(data_train2b)[0]
        data_train2 = numpy.append(data_train2, data_train2b, axis=0)
        
        predicted = train_random_forest(data_train, target_train, data_test, N_trees)
        predicted2 = train_random_forest(data_train2, target_train2, data_test2, N_trees)
        
        #print [x[1] for x in predicted]
        #print [x[1] for x in predicted2]
        predicted3 = combine_predictions(predicted, predicted2)
        #predicted3 = predicted   #one step!
        
        for x in xrange(len(test_index)):
            if target[test_index[x]] == 1:
                predictions_pos[test_index[x]] = predicted3[x][1]
                #print '\t'.join(map(str, ['pos', x, test_index[x], target[test_index[x]], predicted[x,1], predicted2[x,1], predicted3[x][1]]))
            else:
                #print '\t'.join(map(str, ['neg', x, test_index[x], target[test_index[x]], predicted[x,1], predicted2[x,1], predicted3[x][1]]))
                predictions_neg[test_index[x]-len_n] = predicted3[x][1]

    return predictions_pos, predictions_neg



def choose(data, n):
    return data[[numpy.random.choice( data.shape[0], n, replace=False)]]


def getclass_bothpromoters(pos2a, neg2a, pos2b, neg2b):
    n = min(pos2a.shape[0], pos2b.shape[0], neg2a.shape[0], neg2b.shape[0])
    target = [1]*n*2 + [0]*n*2
    
    data1, _ = join_and_balance(choose(pos2a, n), choose(pos2b, n), balance=False)
    data2, _ = join_and_balance(choose(neg2a, n), choose(neg2b, n), balance=False)
    data, _ = join_and_balance(data1, data2, balance=False)
    data = my_transpose(data)[0]
    return get_random_forest(data, target, N_trees)
    


def do_cross_validation(pos1, neg1, pos2, neg2, k, algorithm, estimators, name, outdir):
    #if len(data) < 20*k:
        #k = len(data)/20
    #TODO?
    
    len_n = min(pos1.shape[0], pos2.shape[0], neg1.shape[0], neg2.shape[0])
    target = [1]*len_n + [0]*len_n
    
    data1, data2 = join(pos1, neg1, pos2, neg2)

    kf = cross_validation.StratifiedKFold(target, k)
    auc_list=[]
    mean_tpr = 0.0
    mean_fpr = linspace(0, 1, 100)
    predictions = []
    mean_roc = []   
    cut_value = 0.0

    for train_index, test_index in kf:
        data_train = [data1[index] for index in train_index]
        data_test = [data1[index] for index in test_index]
        
        data_train2 = [data2[index] for index in train_index]
        data_test2 = [data2[index] for index in test_index]
        
        target_train = [target[index] for index in train_index]
        target_test = [target[index] for index in test_index]
        if algorithm == 'svm':
            #TODO
            predicted = train_svm(data_train, target_train, data_test)
        else:
            predicted = train_random_forest(data_train, target_train, data_test, estimators)
            predicted2 = train_random_forest(data_train2, target_train, data_test2, estimators)
            predicted = combine_predictions(predicted, predicted2)
        cut_value += cut_at(target_test, predicted)
        
        roc = pyroc.ROCData([(target_test[i], predicted[i][1],) for i in xrange(0, len(predicted))])
        predictions.append(roc)
        mean_roc += [(target_test[i], predicted[i][1],) for i in xrange(0, len(predicted))]
        auc_list.append(roc.auc())
    mean_tpr /= k
    #mean_tpr[-1] = 1.0
    #mean_auc = metrics.auc(mean_fpr, mean_tpr)
    mean_auc = sum(auc_list)/len(auc_list)
    predictions.append(pyroc.ROCData(mean_roc, 'r-'))
    print("Averaged AUC: %f" % mean_auc)
    print "Averaged cut_value: %f" % (cut_value/k)
    
    
    fig_title = 'ROC Curve for %s on %s \n (mean area = %0.2f)' % (algorithm, name, mean_auc)
    #plot_roc(predictions, fig_title, outdir)
    return mean_auc, cut_value/k


#def save_class(data, target, algorithm, estimators, name, outdir):
    #if algorithm=="random_forest":
        #clf = RandomForestClassifier(n_estimators=estimators)
        #clf.fit(data, target)
        #pickle.dump(clf, open(outdir+name+"_class.pickle", 'w'))     
        #predicted = clf.predict_proba(data)
        #return auc(target, predicted), cut_at(target, predicted)
    

#def plot_roc(predictions, fig_title, outdir):
    #plt.clf()
    #roc_labels=[str(k) for k in xrange(1, len(predictions)+1)]
    #roc_labels[-1]='mean'
    #pyroc.plot_multiple_roc(predictions, title=fig_title, labels=roc_labels, include_baseline=True)
    
    #for i in xrange(1,N_repeats+1):
        #if not os.path.exists(outdir+"roc%d.pdf"%i):
            #plt.savefig(outdir+"roc%d.pdf"%i)
            #break


