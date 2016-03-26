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
#from glob import glob
#from random import sample
from matplotlib import pyplot as plt

import pyroc
import pickle
import sys
import os

from read_data import read_data, join_and_balance
import seaborn as sns


#def train_svm(data, target, test_data, gamma_par=0.01, C_par=32):
    #clf = svm.SVC(kernel='rbf', gamma=gamma_par, C=C_par, probability=True)
    #clf.fit(data, target)
    #return clf.predict_proba(test_data)

def train_random_forest(data, target, test_data, k):
    clf = RandomForestClassifier(n_estimators=k)
    clf.fit(data, target)
    
    predicted = clf.predict_proba(test_data)
    return predicted
    
def get_random_forest(data, target, k):
    clf = RandomForestClassifier(n_estimators=k)
    clf.fit(data, target)
    return clf

def cut_at(target, predicted): #compute cutoff with smallest error    
    npos = len([x for x in target if x == 1])
    nneg = len([x for x in target if x == 0])
    
    #print "CUT AT:"
    #print roc.data, npos, nneg
    
    cut_value = 0.5
    FPR = 0.1
    error = 1
    for cutme in numpy.linspace(1, 0, 21):
        fp = len([1.0 for x in xrange(len(predicted)) if (predicted[x][1] >= cutme) and target[x] == 0 ])
        fn = len([1.0 for x in xrange(len(predicted)) if (predicted[x][1] < cutme) and target[x] == 1 ])
        tmperror = float(fp + fn)/(npos+nneg)
        if tmperror < error:
            error = tmperror
            cut_value = cutme
    fp = len([1.0 for x in xrange(len(predicted)) if (predicted[x][1] >= cut_value) and target[x] == 0 ])
    fn = len([1.0 for x in xrange(len(predicted)) if (predicted[x][1] < cut_value) and target[x] == 1 ])
    print "cutting value %f, error %f, false positives %d, false negatives %d" % ( cut_value, 
                    error, fp, fn )
    return cut_value



def predict_cv(data, k, target):
    kf = cross_validation.StratifiedKFold(target, cv_folds)
    predictions = [None]*len(target)

    for train_index, test_index in kf:
        data_train = data[[train_index]]
        data_test = data[[test_index]]
        target_train = [target[index] for index in train_index]
        target_test = [target[index] for index in test_index]
        predicted = train_random_forest(data_train, target_train, data_test, k)
        
        for x in xrange(len(test_index)):
            predictions[test_index[x]] = predicted[x,1] 
        
    return [predictions[x] for x in xrange(len(target)) if target[x] == 1 ], [predictions[x] for x in xrange(len(target)) if target[x] == 0 ]


def auc(target, predicted):
    roc = pyroc.ROCData([(target[i], predicted[i, 1],) for i in xrange(0, len(predicted))])
    return roc.auc()

def auc2(pos, neg):
    roc = pyroc.ROCData([(1, pos[i],) for i in xrange(0, len(pos))] + 
                            [(0, neg[i],) for i in xrange(0, len(neg))]
                            )
    return roc.auc()

def do_cross_validation(data, k, target, algorithm, estimators, name, outdir):

    kf = cross_validation.StratifiedKFold(target, k)
    auc_list=[]
    mean_tpr = 0.0
    mean_fpr = linspace(0, 1, 100)
    rocs = []
    mean_roc = []   
    cut_value = 0.0

    predictions = [None]*len(target)
    for train_index, test_index in kf:
        data_train = [data[index] for index in train_index]
        data_test = [data[index] for index in test_index]
        target_train = [target[index] for index in train_index]
        target_test = [target[index] for index in test_index]
        if algorithm == 'svm':
            predicted = train_svm(data_train, target_train, data_test)
        else:
            predicted = train_random_forest(data_train, target_train, data_test, estimators)
            
        cut_value += cut_at(target_test, predicted)
        for x in xrange(len(test_index)):
            predictions[test_index[x]] = (test_index[x], predicted[x,1], target_test[x])
        
        roc = pyroc.ROCData([(target_test[i], predicted[i, 1],) for i in xrange(0, len(predicted))])
        rocs.append(roc)
        mean_roc += [(target_test[i], predicted[i, 1],) for i in xrange(0, len(predicted))]
        auc_list.append(roc.auc())
        
        
        
    mean_tpr /= k
    #mean_tpr[-1] = 1.0
    #mean_auc = metrics.auc(mean_fpr, mean_tpr)
    mean_auc = sum(auc_list)/len(auc_list)
    rocs.append(pyroc.ROCData(mean_roc, 'r-'))
    print("Averaged AUC: %f" % mean_auc)
    print "Averaged cut_value: %f" % (cut_value/k)
    
    save_predictions(predictions, outdir)
    
    
    fig_title = 'ROC Curve for %s on %s \n (mean area = %0.2f)' % (algorithm, name, mean_auc)
    plot_roc(rocs, fig_title, outdir)
    
    
    return mean_auc, cut_value/k


def save_class(data, target, algorithm, estimators, name, outdir):
    if algorithm=="random_forest":
        clf = RandomForestClassifier(n_estimators=estimators)
        clf.fit(data, target)
        pickle.dump(clf, open(outdir+name+"_class.pickle", 'w'))     
        predicted = clf.predict_proba(data)
        return auc(target, predicted), cut_at(target, predicted)
    

def save_predictions(predicted, outdir):
    for i in xrange(1,N_repeats+1):
        if not os.path.exists(outdir+"predictions%d.txt"%i):
            with open(outdir+"predictions%d.txt"%i, 'w') as f:
                for i, j, k in predicted:
                    f.write("{0}\t{1}\t{2}\n".format(i,j,k))
            break
            


def plot_roc(predictions, fig_title, outdir):
    plt.clf()
    roc_labels=[str(k) for k in xrange(1, len(predictions)+1)]
    roc_labels[-1]='mean'
    pyroc.plot_multiple_roc(predictions, title=fig_title, labels=roc_labels, include_baseline=True)
    
    for i in xrange(1,N_repeats+1):
        if not os.path.exists(outdir+"roc%d.pdf"%i):
            plt.savefig(outdir+"roc%d.pdf"%i)
            break


