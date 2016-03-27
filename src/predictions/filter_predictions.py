    
import os
from Bio import SeqIO

from ..shared import *
from copy import copy
from ..promoters import promoters
from ..promoters.read_wiggle import *



def save_predictions(regions, filename):
    dirname = os.path.dirname(filename)
    if not os.path.exists(DATAPATH+dirname):
            print 'creating directory', DATAPATH+dirname
            try:
                os.makedirs(DATAPATH+dirname)
            except:
                pass

    with open(DATAPATH+filename, 'w') as f:
        for chrom, region_list in regions.iteritems():
            for region in region_list:
                f.write(str(region)+'\n')


def find_predictions(wiggle, cutoff):
    #returns dictionary with chromosomes and lists of Regions (join regions with high predictions)
    results = {}
    for chrom, regions in wiggle.iteritems():
        result = []
        for region in regions:
            if region.value > cutoff:
                if result != []:
                    last = result[-1]
                    if last.end >= region.start:
                        result[-1].end = region.end
                        continue
                result.append(copy(region))
        results[chrom] = result
    return results

    
def find_best_predictions(wiggle, number, filename):
    #returns dictionary with chromosomes and lists of Regions (join regions with high predictions)
    results = {}
    values = []
    
    count = 0
    
    dirname = os.path.dirname(filename)
    if not os.path.exists(DATAPATH+dirname):
            print 'creating directory', DATAPATH+dirname
            try:
                os.makedirs(DATAPATH+dirname)
            except:
                pass
    
    for chrom, regions in wiggle.iteritems():
        values2=[(region.value, region.chromosome, region.start) for region in regions]
        count += len(values2)
        values += sorted(values2, reverse=True)[:number]
    selected=sorted(values, reverse=True)[:number]
    
    if len(selected) > 1:
    	print "range of predictions of selected windows:", selected[0][0], selected[-1][0]
    print "number of windows (overall):", count
    
    with open(DATAPATH+filename, 'w') as f:
        for chrom, region_list in wiggle.iteritems():
            for region in region_list:
                if (region.value, chrom, region.start) in selected:
                    f.write(str(region)+'\n')
        



if __name__ == '__main__':
    import sys
    print sys.argv
    if len(sys.argv) <= 2:
        print "USAGE: predict.py input_filename output_name threshold|best [threshold|N]"
        sys.exit(1)
    
    wiggle=read_all(sys.argv[1])
    
    if sys.argv[3] == 'threshold':
        threshold=float(sys.argv[4])
        save_predictions(find_predictions(wiggle, threshold), sys.argv[2])
    if sys.argv[3] == 'best':
        find_best_predictions(wiggle, int(sys.argv[4]), sys.argv[2])

        
    
