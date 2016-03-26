
from ..shared import *
import load_vista

import sys








if  __name__ =='__main__':
    
    if len(sys.argv) < 4:
        print "USAGE: add_sequences.py TYPE inputfile outputfile [start_id] [tissue]"
        print "Genome file: %s%s, files_directory: %s"% (DATAPATH, GENOME_FILE, DATAPATH)
        sys.exit(1)
    
    if len(sys.argv) == 6:
        tissue = [sys.argv[5]]
        positive = True
    else:
        tissue = []
        positive = False
    
    if sys.argv[1] == 'bed':
        positions = load_vista.load_bed(sys.argv[2], positive, tissue, startid=int(sys.argv[4]))
        load_vista.add_seqs(positions, sys.argv[3])
    else: #vista type
        positions = load_vista.load_enhancers(sys.argv[2])
        load_vista.add_seqs(positions, sys.argv[3])
        