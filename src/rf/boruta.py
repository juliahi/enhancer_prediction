



from ..shared import *

import pickle 

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import DataFrame, FloatVector, IntVector
#from rpy2.robjects.numpy2ri import numpy2ri
import rpy2.robjects.numpy2ri


def start_boruta():
    #konfiguracja rpy potrzebna przed uruchomieniem Boruty. Nie moze byc wykonana wielokrotnie
    robjects.conversion.py2ri = rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()


def get_decisions(result):
    boruta = importr('Boruta')
    r = robjects.r
    
    stats = boruta.attStats(result)
    return [(name, float(meanZ), stats.rx2(6).levels[dec-1]) for (name, meanZ, dec) in zip(stats.rownames, stats.rx2(1), stats.rx2(6))]

def run_boruta(data, target, names, name, outdir):
    #uruchomienie algorytmu Boruta na data i target
    grdevices = importr('grDevices')
    boruta = importr('Boruta')
    r = robjects.r
    base = importr('base')

    data2 = {}
    for i in xrange(len(names)):
        data2[names[i]] = FloatVector((data[:,i]))

    x = robjects.DataFrame(data2)
    y = IntVector((target))
    print "running Boruta"
    result = boruta.Boruta(x, y)
    print result
    print boruta.attStats(result)
    f = file(outdir+"boruta.data", 'w')
    pickle.dump(result, f)
    f.close()
    

    #grdevices.pdf(outdir+'boruta.pdf', width=60, height=30)
    #graphics = importr('graphics')
    #graphics.par(**{'cex.axis': 1, 'cex.lab': 2, 'mar': [16, 8, 4, 2], 'mgp': [14, 1, 0]})
    #r.plot(result, ylab = '', las=2)
    #graphics.mtext('Importance', side=2, line=5, cex=2)
    #grdevices.dev_off()
    #r['dev.off']
    
