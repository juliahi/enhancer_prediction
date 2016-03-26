"""
An implemmentation of a method for comparing AUCs of different classifiers based on the same dataset.

method from JAmes A. Hanley and Barbara J. McNeil, published in Radiology, Vol 148, No 3, pp 839-843, 1983.
"""

import math
from rpy2 import robjects


def seROC(AUC,pos,total):
    """
    Estimate of standard error for a classifier with a given AUC

    parameters:
    AUC = a value of AUC (0.0-1.0)
    pos - number of positive cases in the sample
    total - number of all cases in the sample
    """
    na=pos
    nn=total-pos
    a=AUC
    a2=a**2
    q1=a/(2-a)
    q2=(2*a2)/(1+a)
    se=math.sqrt((a*(1-a)+(na-1)*(q1-a2)+(nn-1)*(q2-a2))/(nn*na))
    return se

def zp(z):
    """p-value (two tailed) for a z-score"""
    return robjects.r.dnorm(z)[0]

def pAUC(AUC,pos,total):
    """
    p-value of an AUC score against H0="classifier is random"
    
    AUC = a value of AUC (0.0-1.0)
    pos - number of positive cases in the sample
    total - number of all cases in the sample
    """
    se=seROC(0.5,pos,total)
    return zp((AUC-0.5)/se)

def diffAUC(AUC_TEST,AUC_EXPECT,pos,total):
    """
    p-value of an AUC score against H0="classifier AUC_TEST is no better than AUC_EXPECT"
    
    AUC_TEST,AUC_EXPECT = a value of AUC (0.0-1.0)
    pos - number of positive cases in the sample
    total - number of all cases in the sample
    """
    se=seROC(AUC_EXPECT,pos,total)
    return zp((AUC_TEST-AUC_EXPECT)/se)


def cROC(AUC1,AUC2,pos,total,r):
    """
    p-value of the difference between two classifiers given their AUCs and their correlation

    H0=there is no difference between the strength of classifiers giving AUC1 and AUC2

    AUC1,AUC2 = a value of AUC (0.0-1.0)
    
    pos - number of positive cases in the sample
    total - number of all cases in the sample

    r - r-coefficient of pearson correlation for vectors of predictions of the two classifiers
    """
    
    se1=seROC(AUC1,pos,total)
    se2=seROC(AUC2,pos,total)
    sed=math.sqrt(se1**2+se2**2-2*r*se1*se2)
    return zp((AUC1-AUC2)/sed)


###EXAMPLE OF USAGE FOLLOWS

def read_pred(f):
    import collections
    res=collections.defaultdict(dict)

    for ln in open(f):
        if ln[0]==">":
            tis=ln.strip()[1:].split("\t")[0]
            pred=res[tis]
        else:
            rec=ln.strip().split("\t")
            if rec[2]=="True":
                lab=1
            else:
                lab=0
            pred[rec[0]]=(float(rec[1]),lab)
                
    return res


if __name__=="__main__":
    import sys,Rocr
    from numpy import array
    from rpy2 import robjects
    import rpy2.robjects.numpy2ri

    p1=read_pred(sys.argv[1])
    p2=read_pred(sys.argv[2])
    for tis in p1.keys():
        keys=[g for g in p1[tis] if g in p2[tis]]
        #print p1.keys(),p1[tis].keys()[:10]
        pred1=[p1[tis][k][0] for k in keys]
        lab1=[p1[tis][k][1] for k in keys]
        pred2=[p2[tis][k][0] for k in keys]
        lab2=[p2[tis][k][1] for k in keys]
        #print len(pred1),len(pred2)
        r=robjects.r.cor(array(pred1),array(pred2))[0]
        pos=len([1 for k in keys if p1[tis][k][1]==1])
        total=len(pred1)
        a1=Rocr.aucROC(pred1,lab1)
        a2=Rocr.aucROC(pred2,lab2)
        
        print tis,a1,pAUC(a1,pos,total),a2,pAUC(a2,pos,total),r,cROC(a1,a2,pos,total,r)
