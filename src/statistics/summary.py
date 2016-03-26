



from ..shared import *

import roc_compare
import pandas

import itertools
import numpy


#DB       POS     NEG     kmers   hmods   dist    ntrees  use_GC  cv_folds        used_class      AUC     date    cut_at  output_dir      POS_size        NEG_size

def filtr(df):
    return df[(df.dist == True) & (df.ntrees ==100) & (df.use_GC == False)].drop(['output_dir', \
            'date', 'dist', 'ntrees', 'use_GC', 'used_class'], axis=1)
#DB       POS     NEG     kmers   hmods   cv_folds       AUC     cut_at    POS_size        NEG_size


def add_erwin(results):
    
    n = len(results)
    for tissue, kmers, hmods, auc, npos in [ ['heart',"['4mers']", 'ALL', 0.85, 84],
        ['heart', '-', 'ALL', 0.86, 84],
        ['heart',"['4mers']", 'No', 0.83, 84],
        
        #['brain',"['4mers']", 'ALL',  0.73, ? ],
        
        ['both',"['4mers']", 'ALL',  0.96, 711],
        ['both', '-', 'ALL',  0.89, 711],
        ['both',"['4mers']", 'No', 0.88, 711],
]:
        results.loc[n] = ['erwin', tissue, 'randoms_'+tissue, kmers, hmods, 10, auc, None, npos, npos ]
        #add ERWINs results
        n += 1
    
    

def read_summary(db=None, tissue=None, erwin=True):
    
    
    summary = pandas.read_csv(SUMMARYFILE, sep='\t')
    summary = filtr(summary)
    if erwin:
        add_erwin(summary)
    
    groups = summary.groupby(('DB', 'POS', 'kmers', 'hmods')).agg({'AUC': [numpy.mean, numpy.std, len], 'POS_size': max, 'NEG_size': max})
    
    if db:
        if db in list(groups.index.levels[0]):
            groups2 = groups.loc[db]
        else:
            groups2 = groups.loc["['"+db+"', '"+db+"']"]
        if db == 'vista' and erwin:
            groups = groups2.append(groups.loc['erwin'])
        else:
            groups=groups2
    if tissue:
        groups = groups.loc[tissue]
        
    
        
    return groups
        



if __name__ == '__main__':
    

    s = read_summary('vista', erwin=False)
    print s
    
    dhmods = dict(zip(['-', 'h1hesc.list', 'tier12.list'], ['-', 'H1hESC', 'Tier1\&2']))
    
    f = open('results_kmers.csv', 'w')
    f.write("tissue,kmers,hmods,AUC\n")
    for kmers in ["['3mers']", "['4mers']"]:
        for tissue in ['heart', 'brain']:
            #print s.loc[tissue, kmers]['hmods']
            for hmods in ['-', 'h1hesc.list', 'tier12.list']:
                try:
                    #print tissue, kmers, hmods, s.loc[tissue, kmers, hmods]['AUC']['mean']
                    val=s.loc[tissue, kmers, hmods]['AUC']['mean']
                    f.write("{0},{1},{2},{3:.3f}\n".format(tissue, kmers.strip("[]'"), dhmods[hmods], val))
                except:
                    pass
    f.close()
    
    f = open('results_tissues.csv', 'w')
    f.write("tissue,kmers,hmods,AUC\n")
    for tissue in ['both', 'heart', 'brain']:
        for kmers in ["['3mers']", "['4mers']"]:
            #print s.loc[tissue, kmers]['hmods']
            for hmods in ['-', 'h1hesc.list', 'tier12.list']:
                try:
                    #print tissue, kmers, hmods, s.loc[tissue, kmers, hmods]['AUC']['mean']
                    val=s.loc[tissue, kmers, hmods]['AUC']['mean']
                    if tissue == 'both':
                        f.write("non-specific,{0},{1},{2:.3f}\n".format( kmers.strip("[]'"), dhmods[hmods], val))
                    else:
                        f.write("{0},{1},{2},{3:.3f}\n".format(tissue, kmers.strip("[]'"), dhmods[hmods], val))
                except:
                    pass
    f.close()