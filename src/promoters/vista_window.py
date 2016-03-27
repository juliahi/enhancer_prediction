from Bio import SeqIO

from ..rf import vista
from ..shared import *
from ..preprocessing import load_vista
import promoters

import sys

def gen_windows(length, outname, db):
    if db == 'vista':
        enhancers = vista.load_target("all", False)
    elif db == 'randoms':
        enhancers = vista.load_target("randoms_heart", False) + vista.load_target("randoms_brain", False)
    else:
        enhancers = vista.load_target(db, False) 
        
    
    f = open(DATAPATH+outname, 'w')
    record_dict = SeqIO.index(DATAPATH+GENOME_FILE, "fasta")
    print "genome loaded"
    for e in enhancers:
        olds =  e.start 
        olde = e.end
        e.start = max(0, e.start + (e.end-e.start)/2 - length/2)
        e.end = e.start + length - 1        #endpoint included
        if olds > e.start or olde < e.end:
            e.seq = str(record_dict[e.chromosome].seq[e.start:(e.end+1)])
        else:
            e.seq = e.seq[(e.start-olds):(e.end-olde)]
        assert len(e.seq)==length, "length of seq not correct! %d != %d"%(len(e.seq), length)
        f.write(str(e)+"\n")

if __name__ == '__main__':    
    
    if len(sys.argv) <= 3:
        print "USAGE: vista_window.py window_size outname db"
        sys.exit(1)
        
    window_size = int(sys.argv[1])
    outname = sys.argv[2]
    db = sys.argv[3]
    gen_windows(window_size, outname, db)
