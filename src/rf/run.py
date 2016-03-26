# -*- coding: utf-8 -*-

#import matplotlib
#matplotlib.use('Agg')


import boruta
import vista
import fantom
import new_train as train
from ..shared import *
from ..promoters import promoters

import argparse
import sys
import os
import random
import pickle
import numpy
import pyroc
import time

from matplotlib import pyplot as plt

from load import load_data, load_target
from read_data import cut_rows, my_transpose

from glob import glob
    



def predict(datapos, dataneg, class_dir, outdir):
    
    try:
        class_filename = glob(RESULTSPATH+class_dir+'/*')[0]
    except:
        assert True, 'No classifier in %s'%class_dir
        
    print "Predicting using classifier from ", class_filename
    classifier = pickle.load(open(class_filename))
    
    data, target = train.join_and_balance(datapos, dataneg, False)
    data, names = my_transpose(data)
    predicted = classifier.predict_proba(data)
    roc = pyroc.ROCData([(target[i], predicted[i, 1],) for i in xrange(0, len(predicted))])
    
    train.save_predictions([(i, predicted[i, 1],target[i]) for i in xrange(0, len(predicted))], outdir)
    
    print "AUC=", roc.auc()
    #plot roc ?
    return [roc.auc(), 0]
    
    
def get_class(datapos, dataneg):
    data, target = train.join_and_balance(datapos, dataneg)
    data, names = my_transpose(data)
    clas = train.get_random_forest(data, target, N_trees)
    return clas
    
def get_class2(data, target):
    clas = train.get_random_forest(data, target, N_trees)
    return clas
    
def train_cv(datapos, dataneg, name, outdir, boruta_arg):
    len_neg = len(dataneg)
    auc = 0
    cut_at = 0.0
    
    for x in xrange(N_repeats):
        if len_neg > len(datapos)*(N_repeats-1):
            print "splitting negatives into %d groups" % N_repeats
            pocz = int(round(x*len_neg/N_repeats))
            kon = int(round((x+1)*len_neg/N_repeats))
            data, target = train.join_and_balance(datapos, cut_rows(dataneg, range(pocz,kon)))
        else:
            data, target = train.join_and_balance(datapos, dataneg)
        data, names = my_transpose(data)
        result = train.do_cross_validation(data, cv_folds, target, "random_forest", N_trees, name, outdir)
        auc += result[0]
        cut_at += result[1]
    
    
    
    print "Mean AUC in %d repeats = %f" % (N_repeats, auc/N_repeats)
    print "Mean cut value = %f" % (cut_at/N_repeats)
    
    
    
    sys.stdout.flush()
    if boruta_arg:
        boruta.run_boruta(data, target, names, name, outdir)
    return auc/N_repeats, cut_at/N_repeats
        
def train_save(datapos, dataneg, name, outdir):
    data, target = train.join_and_balance(datapos, dataneg)
    data, names = my_transpose(data)
    
    auc = train.save_class(data, target, "random_forest", N_trees, name, outdir)
    sys.stdout.flush()
    #if boruta_arg:
    #    boruta.run_boruta(data, target, names, name, outdir)
    #return auc
        

tissue_choices = ('heart', 'brain', 'both', 'positives', 'balanced',
                      'randoms_heart', 'randoms_brain', 'randoms_both', 'randoms_positives', 
                      'predicted_promoters_4mers_heart_0.8', 
                      'predicted_promoters_h1hesc_4mers_heart_0.8',
                      'predicted_promoters_4mers_brain_0.8', 
                      'predicted_promoters_h1hesc_4mers_brain_0.8',
                      'predicted_promoters_4mers_both_0.8',
                      'predicted_promoters_h1hesc_4mers_both_0.8',
                      'predicted_promoters_4mers_balanced_0.8',
                      'predicted_promoters_h1hesc_4mers_balanced_0.8',
                      'heart_notss', 'brain_notss', 'both_notss', 'positives_notss',
                      
                      'predicted_promoters_4mers_both_0.0',
                      )

def main(argv):
    
    parser = argparse.ArgumentParser(description='Run test')
    parser.add_argument('-o', '--output', type=str, 
                   help='output directory')
    
    parser.add_argument('--db', nargs='+',  choices=('vista', 'vista1500', 'fantom', 'fantom_long', 'promoters'), required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--saveclass', action='store_true', help='train classifier on whole dataset and save to pickle')
    group.add_argument('--useclass', type=str, help='use classifier saved in file')
    
    parser.add_argument('-p', '--positives', choices=tissue_choices, required=True)
    parser.add_argument('-n', '--negatives',  choices=tissue_choices, required=True)
    parser.add_argument('--boruta', action='store_true')
    parser.add_argument('--distinct', action='store_true')
    parser.add_argument('--histmods', type=str, default='', help='histone modificiations list file')
    parser.add_argument('--kmers', type=str, nargs='+', default='', help='kmers extension')
    
    #parser.add_argument('--usegc', action='store_true')
    
    args = parser.parse_args()
    
    if args.kmers == None and args.histmods == None:
        parser.error('Specify kmers or histmods')
    
    if len(args.db) == 1:
        args.db.append(args.db[0])
    
    if 'fantom' in args.db and ('positives' in args.positives or 'positives' in args.negatives ):
        parser.error('Positives for FANTOM not defined')
    
    
    #prepare output directory
    outdir = RESULTSPATH+args.output

    try:
        
            maks = 0
            for i in xrange(100, 0, -1):
                if os.path.exists(outdir+'_'+str(i)):
                    maks = i
                    break
            #print 'moving to', outdir+str(maks+1)
            
            outdir += '_'+str(maks+1)
            os.mkdir(outdir)
    except:
        parser.error('Cannot create directory %s' % outdir)
    outdir += '/'
    
    #write report, redirect stdout to log file
    orig_stdout = sys.stdout
    outfile = open(outdir+'log.txt', 'w')
    sys.stdout = outfile
    
    print args
    #print "cv_folds=%d, N_trees=%d, N_repeats=%d" % ( shared.cv_folds, shared.N_trees, shared.N_repeats )
    print "cv_folds=%d, N_trees=%d, N_repeats=%d, usegc=%s" % ( cv_folds, N_trees, N_repeats, USEGC )

    
    #load pos and neg data without balance
    datan = load_data(args.db[1], args.histmods, args.kmers, args.negatives, args.distinct )
    datap = load_data(args.db[0], args.histmods, args.kmers, args.positives, args.distinct )
    print "data sizes: %s %d %d, %s %d %d" % (args.positives, datap.shape[0], len(datap.dtype.names), args.negatives, datan.shape[0], len(datan.dtype.names))
    sys.stdout.flush()
    if args.boruta:
        boruta.start_boruta()
        
    
    #load classifier from pickle
    if args.useclass:
        auc = predict(datap, datan, args.useclass, outdir )
    #train new classifier
    else:
        if args.saveclass:
            name = args.db[0] + '-' + args.db[1]+"_" + args.positives + "_vs_" + args.negatives + "_" + args.histmods + "_" + str(args.kmers)
            auc = train_save(datap, datan, name, outdir)
            
        #else:
        name = args.positives + " vs " + args.negatives
        auc = train_cv(datap, datan, name, outdir, args.boruta)
        
    
    #finish 
    #all options: DB, POS, NEG, kmers, hmods, ntrees, cv_folds (0 if no cv), used_class, auc
    
    if args.kmers == '': args.kmers = '-'
    if args.histmods == '': args.histmods = '-'
    summary = "%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s" % (args.db, args.positives, 
                                              args.negatives, args.kmers, 
                                              args.histmods, args.distinct, N_trees, USEGC)
    
    if not args.useclass and not args.saveclass:
        summary += "\t%d\t" % (cv_folds )
    elif args.useclass:
        summary += "\t0\t%s" % (args.useclass)
    else:
        summary += "\t0\t%s" % (args.saveclass)
    summary += "\t%f" % auc[0]
    summary += "\t%s" % time.ctime()
    summary += "\t%f" % auc[1]
    summary += "\t%s" % outdir
    summary += "\t%d\t%d" % (datap.shape[0], datan.shape[0])
    
    summary += "\n"
    
    print summary
    
    sys.stdout = orig_stdout
    outfile.close()
    
    summaryf = open(SUMMARYFILE, 'a+')
    summaryf.write(summary)
    summaryf.close()
    

if __name__ == "__main__":
   main(sys.argv[1:])

