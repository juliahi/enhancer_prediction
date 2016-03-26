#check if there are overlaps btw groups of sequences


from load_vista import load_enhancers, filter_type
import sys

from ..shared import *
from ..rf.load import load_target

def overlap(data, chromosome, start, end):
    for x in data:
        if x.chromosome == chromosome:
            if ((x.start <= end) and (start <= x.end)):
                return True
    return False



def check_overlap(data1, data2):
    
    d = {}
    for x in data2:
        if d.has_key(x.chromosome):
            d[x.chromosome].append((x.start, x.end, x))
        else:
            d[x.chromosome] = [(x.start, x.end, x)]
    overlaps = []
    for x in data1:
        if d.has_key(x.chromosome):
            for (start, end, info) in d[x.chromosome]:
                if ((x.start <= end) and (start <= x.end)):
                    print ("%d and %d: %s: %d - %d, %d - %d" % (x.id, info.id, x.chromosome, start, end, x.start, x.end))
                    overlaps.append(x.id)
                    break
                
    overlaps=list(set(overlaps))
    return overlaps
 
 
 

if __name__ == "__main__": 
    dist = True
    if len(sys.argv) < 4:
        print "USAGE: overlaps.py from_file with_db tissue output"
        print "output_directory: %s"% ( DATAPATH)
        sys.exit(1)
    
    file1 = sys.argv[1]
    db = sys.argv[2]
    tissue = sys.argv[3]
    outname = sys.argv[4]
    
    
    remove = sorted(check_overlap(load_enhancers(file1), load_target(db, tissue, False)))
    print remove
    
    with open(DATAPATH+outname, 'a+') as f:
        
        #remove = [x  for x in remove]
        #print "removing %d"% len(set(remove))
        f.write("\n".join([str(x) for x in remove])+"\n")
    f.close()

    
