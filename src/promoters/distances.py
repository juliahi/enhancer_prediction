

from ..shared import *
from promoters import *
from ..rf import run
from ..rf import load
import sys
import read_wiggle
import predict
import pickle

def dist(seq, tss, direction):
    
    if seq[0] <= tss <= seq[1]:
        return 0, tss
    if seq[1] < tss:
        return -direction*(tss-seq[1]), tss
    if seq[0] > tss:
        return direction*(seq[0]-tss), tss


def bins(val, window):
    if val == 0: return 0
    if 0 < val <= window: return 1
    if 0 > val >= -window: return -1
    if val > window: return 2
    if val < -window: return -2



def count_dists(sequences, promoters, dist_fun): 
    if sequences == []:
        return []
    if type(sequences) == list :
        return count_dists_lists(sequences, promoters, dist_fun)
    elif type(sequences) == dict:
        return count_dists_regions(sequences, promoters, dist_fun)
    
def count_dists_lists(sequences, promoters, dist_fun): 
    dists = []
    
    for s in sequences:
            d,pos = distance(s.start, s.end, promoters[s.chromosome], dist_fun)
            #if write overlap with promoters in BED format:
            #if d == 0:
                #print "%s\t%d\t%d\tid=%d,TSS=%d\t0\t."% (s.chromosome, s.start, s.end, s.id, pos)
            dists.append(bins(d,WINDOW))
    print "-1: %d, 0: %d, 1: %d" %(dists.count(-1), dists.count(0), dists.count(1))
    print "%d from %d sequences in promoters of length %d"%(len([ 1 for x in dists if -1 <= x <=0]), len(dists), WINDOW)
    return dists


#def distance(dists):
    #dists1 = [x for x in dists if x >= 0]
    #if dists1 != []:
        #tmp = min( dists1 )
        #if tmp == 0: return 0
    #dists2 = [x for x in dists if x < 0]
    #if dists2 != []:
        #tmp2 = -max( dists2 )
    #else:
        #return tmp
    #if tmp < tmp2:
        #return tmp
    #else:
        #return -tmp2

def distance(start, end, promoters, dist_fun):
    dminus = None
    dplus = None
    for x in promoters:
            d = dist_fun((start, end), x[0], x[1]) 
            d = d[0]
            if d == 0:
                #dplus = 0
                return 0, x[0]
                #break
            if d > 0:
                if dplus == None or d < dplus:
                    dplus = d
                    ppos = x[0]
            else:
                if dminus == None or d > dminus:
                    dminus = d
                    mpos = x[0]
            #if abs(d) < 5000:
                #print d, start, end, x[0], x[1]
    if dplus == None:
        return dminus,mpos
    if dminus == None:
        return dplus,ppos
    if -dminus < dplus:
        return dminus,mpos
    else:
        return dplus,ppos

def count_dists_regions(sequences, promoters, dist_fun): 
    dists = []
    
    for chrom, regions in sequences.iteritems():
        for region in regions:
            #chrom = region.chromosome
            start = region.start
            end = region.end
  
            d,_ = distance(start, end, promoters[chrom], dist_fun)
            dists.append(bins(d,WINDOW))
            
    print "-1: %d, 0: %d, 1: %d" %(dists.count(-1), dists.count(0), dists.count(1))
    print "%d from %d sequences in promoters of length %d"%(len([ 1 for x in dists if -1 <= x <=0]), len(dists), WINDOW)
    return dists



def count_distances(filename):
    
    tss = read_tss()
    print "Number of promoters:", len(tss)
    
    #vista_heart = load.load_target('vista', 'heart', True) 
    #vista_brain = load.load_target('vista', 'brain', True) 
    #randoms_heart = load.load_target('vista', 'randoms_heart') 
    #randoms_brain = load.load_target('vista', 'randoms_brain') 
    #fantom_heart = load.load_target('fantom', 'heart', True) 

    #erwin_step1 = read_wiggle.read_bed(DATAPATH+'enhancerfinder_step1_hg19.bed') 
    #erwin_heart = read_wiggle.read_bed(DATAPATH+'enhancerfinder_heart_hg19.bed') 
    #erwin_brain = read_wiggle.read_bed(DATAPATH+'enhancerfinder_brain_hg19.bed') 
    
    cutoff = 0.8
    #michal1 = predict.find_predictions(read_wiggle.read_wiggle(WIGGLEPATH+'vista_both__4mers/chr21.id.wig'), cutoff)    
    #michal2 = predict.find_predictions(read_wiggle.read_wiggle(WIGGLEPATH+'vista_both_h1hesc.list_4mers/chr21.id.wig'), cutoff)

    f = open(DATAPATH+filename, 'a+')
    for name, sequences in [
                        #('vista_heart', vista_heart),('vista_brain', vista_brain),
                      #('randoms_heart', randoms_heart),
                      #('randoms_brain', randoms_brain),
                      #('erwin_heart', erwin_heart), 
                      #('erwin_brain', erwin_brain), 
                      #('erwin_step1', erwin_step1),
                      ('predictions vista_4mers_heart > 0.8', 
                                load.load_target('predictions', 'predicted_vista_heart__4mers_0.8')),
                      ('predictions vista_h1hesc_4mers_heart > 0.8', 
                                load.load_target('predictions', 'predicted_vista_heart_h1hesc.list_4mers_0.8')),
                      ('predictions vista_4mers_brain > 0.8',
                                load.load_target('predictions', 'predicted_vista_brain__4mers_0.8')),
                      ('predictions vista_h1hesc_4mers_brain > 0.8', 
                                load.load_target('predictions', 'predicted_vista_brain_h1hesc.list_4mers_0.8')),
                      #('predictions vista_4mers_both > 0.8', 
                                #load.load_target('predictions', 'vista_heart__4mers_0.8'),
                      #('predictions vista_h1hesc_4mers_both > 0.8', 
                                #load.load_target('predictions', 'vista_heart__4mers_0.8'),
                      ('vista_heart1500', load.load_target('vista1500', 'heart', True)), 
                      ('vista_brain1500', load.load_target('vista1500', 'brain', True)), 
                      ('randoms_heart1500', load.load_target('vista1500', 'randoms_heart')),
                      ('randoms_brain1500', load.load_target('vista1500', 'randoms_brain')),  
                       ]:
        print 'Counting for', name, len(sequences)
        dists = count_dists(sequences, tss, dist)
        pickle.dump([name, dists], f)
    f.close()
    #return distances

def load_distances(filename):
    f = open(DATAPATH+filename)
    result = []
    while True:
        try:
            result.append(pickle.load(f))
        except EOFError:
            break
    return result





if __name__ == '__main__':
    
    if len(sys.argv) <= 1:
        print "USAGE: distances.py outfile"
        sys.exit(1)
    
    count_distances(sys.argv[1])
    
    #tss = read_tss()
    #sequences = load.load_target('vista', 'heart', dist=False)
    #count_dists(sequences, tss, dist)