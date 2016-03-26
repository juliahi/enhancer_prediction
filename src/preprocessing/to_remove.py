


import sys
from ..shared import *


def load_fails(name):
    fails_file = open(DATAPATH+name)
    names = fails_file.readline().split('\t')[1:]
    
    n = len(names)
    fails = [0]*n
    n_lines = 0
    for line in fails_file:
            fails = [fails[i] + int(line.split('\t')[i+1]) for i in range(n)]
            n_lines += 1
    fails = [float(x)/n_lines for x in fails]
    return zip(names, fails)


def fails_to_remove(fails_file, out_file):
    # change enhancers with more than 25% fails
    f = open(DATAPATH+out_file, "a+")
    
    for k,v in load_fails(fails_file):
        print k,v
        if v > 0.25:
            f.write("%s\n"%k)
    f.close()




if __name__=="__main__":
    
    if len(sys.argv) < 3:
        print "USAGE: to_remove.py fails_filename removes_filename"
        print "data directory: %s"% (DATAPATH)
        sys.exit(1)
    
    fails_to_remove(sys.argv[1], sys.argv[2])