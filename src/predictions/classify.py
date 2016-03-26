

import predictions
from ..rf.load import load_data
import pickle
import sys
import os

from ..shared import *
from ..rf.read_data import my_transpose

from glob import glob

import sklearn
from sklearn import svm, grid_search
from sklearn import cross_validation
from sklearn import metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn import tree 



def classify(name, class_file, db, tissue, histmods, kmers):
    print sklearn.__version__
    chrom = name.split('/')[2]
    step = name.split('/')[1]
    
    
    if histmods == '-': histmods = ''
    else: histmods += '.list'
    
    data = load_data('whole', histmods, [kmers], chrom)    
    data, names = my_transpose(data)
    classifier = pickle.load(open(class_file, 'r'))
    predicted = classifier.predict_proba(data)
    
    
    
    path=DATAPATH+'predictions/%s/%s_%s_%s_%s' % (step, db, tissue, histmods, kmers)
    #second step
    #path=DATAPATH+'predictions/%s/second_step_%s_%s_%s_%s' % (step, db, tissue, histmods, kmers)
    try:
            results = open(path+'/%s.wig'%chrom, 'w+')
    except:
        try:
            os.makedirs(path)
        except:
            pass
        results = open(path+'/%s.wig'%chrom, 'w+')
    
    trackname = "db=%s tissue=%s histmods=%s kmers=%s chrom=%s"%(db, tissue, histmods, kmers, chrom)
    results.write(('track type=wiggle_0 name="%s" description="enhancers prediction"  visibility=full autoScale=off ')%trackname)
    results.write(('vieLimits=0.0:1.0 color=50,150,255 yLineMark=11.76 yLineOnOff=on priority=10\nvariableStep chrom=%s span=%s\n') % (chrom, step))
    list_adnotation = predictions.load_positions(name)
    
    for pos, pred in zip(list_adnotation, predicted):
        results.write(('%s %f\n') % (pos, pred[1]))

    results.close()
    
    

if __name__ == "__main__":
    if len(sys.argv) <= 5:
        print "USAGE: classify.py predictions/WINDOW/chrN.id db tissue histmods kmers "
        sys.exit(1)
        
    name = sys.argv[1]
    db = sys.argv[2]
    tissue = sys.argv[3]
    histmods = sys.argv[4]
    kmers = sys.argv[5]
    
    
    class_path = 'classifiers/class_%s_%s'%(db, tissue)
    
    if histmods != '-': class_path += '_' + histmods
    if kmers != '-': class_path += '_' + kmers
    
    #second step
    #classpath = 'classifiers/class_balancedpromoters_randoms_both_4mers'
    print class_path
    class_file = glob(RESULTSPATH+class_path+'*/*_class.pickle')
    print class_file
    classify(name, class_file[0], db, tissue, histmods, kmers)
    
    
