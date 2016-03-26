from load_vista import load_enhancers
from ..shared import *
from bx.bbi.bigwig_file import BigWigFile
from glob import glob
import os
import StringIO
import time
import fileinput
import sys
import random

#from multiprocessing import Pool
#from multiprocessing.dummy import Pool as ThreadPool 


def count_chrom_mean(bigwig):
    mean_all = {}
    for (chrom, length) in chr_lens.items():
        if chrom != 'chrY':
            summary = bigwig.summarize(chrom, 1, length, 1)
        
            mean_all[chrom] = summary.sum_data / summary.valid_count 
    return mean_all

def count_mean_signal(enhancers, bigwig_file, name):
    print "processing", bigwig_file



    if os.path.exists(bigwig_file + ".means." + name.split('/')[-1]):
        print "file exists:", bigwig_file+ ".means." + name.split('/')[-1]
        return []
    
    #output = StringIO.StringIO()
    
    
    f = open(bigwig_file, "r") 
    print "bigwig file opened"
    bigwig = BigWigFile(file=f)
    mean_all = count_chrom_mean(bigwig)
    print "chromosome means counted"


    output2 = open(bigwig_file + ".means." + name.split('/')[-1], "w")
    print "output file opened", bigwig_file+".means."+name.split('/')[-1] 
    

    start = time.clock()
    fails = []
    i = 0
    for enh in enhancers:
        if i%10000 == 0:
            print bigwig_file, name, i

        summary = bigwig.summarize(enh.chromosome, enh.start, enh.end+1, 1)   

        #+1 added 24.09.15 after finding endpoint not included
        if summary.valid_count*10 < enh.end-enh.start+1:
            mean = mean_all[enh.chromosome]
            #fails.write("%d\t%d\n" % (enh.id, summary.valid_count))
            fails.append(1)
        else:
            mean = summary.sum_data / summary.valid_count
            fails.append(0)
        i += 1
        output2.write("%d\t%f\n" % (enh.id, mean))
    #bigwig.close()
    #output2 = open(bigwig_file + ".means." + name, "w")
    #output2.write(output.getvalue())
    #output.close()
    output2.close()
    f.close()

    print "output written to: %s.means.%s" % (bigwig_file, name.split('/')[-1])
    end = time.clock()
    print "time: %.2f s" % (end-start)
    return fails
    
    
def count_all_means(enhancers_file, directories, name):
    enhancers = load_enhancers(enhancers_file)
    print "nr of sequences:", len(enhancers)
    names = [str(x.id) for x in enhancers]
    try:
        fails_file = open(DATAPATH+"%s.fails2" %(name), "a")
    except:
        fails_file = open(DATAPATH+"%s.fails2" %(name), "w+")
        fails_file.write('\t'.join([' ']+names)+"\n")
    
    print "saving fails number to %s%s.fails"%(DATAPATH, name)
    n = len(names)
    result = [0]*n
    
    #n_lines = 0
    filenames = []
    N_pools = 4
    for directory in directories:
        filenames += glob(directory + "/*.bigWig")+glob(directory + "/*.bw")
    
    #pool = ThreadPool(N_pools) 
    #fail_list = pool.map(lambda x: count_mean_signal(enhancers, x, name), filenames)
    #pool.close()
    
    random.shuffle(filenames)
    for filename in filenames:
            fails = count_mean_signal(enhancers, filename, name)
            if fails != []:
                fails_file.write("\t".join([filename.split('/')[-1]]+map(str, fails))+"\n")
            #result = [result[i] + fails[i] for i in range(n)]
            #n_lines += 1
    #result = [float(x)/n_lines for x in result]
    #return zip(names, result)


def check_position(chrom, start, end):
    #is there 10% coverage of region [start, end]
    valids = 0.
    wrong = 0.
    for directory in [x[0] for x in os.walk(DATAPATH+"data")]:
        for filename in glob(directory + "/*.bigWig")+glob(directory + "/*.bw"):
            f = open(filename, "r") 
            bigwig = BigWigFile(file=f)
            summary = bigwig.summarize(chrom, start, end+1, 1)   
            if summary.valid_count*10 < end-start+1:
                wrong += 1
            else:
                valids += 1
    return (valids/(valids+wrong) >= 0.75)







if __name__=="__main__":
    
    if len(sys.argv) < 2:
        print "USAGE: histmods.py sequences_file [cell_type]"
        print "data directory: %s"% (DATAPATH)
        sys.exit(1)
    if len(sys.argv) == 3:
        cell_types = [sys.argv[2]]
   
    if cell_types[0] == "X":
        
        count_all_means(sys.argv[1], ['/home/julia/data/'+x for x in cell_types], sys.argv[1])
    else:
        count_all_means(sys.argv[1], [DATAPATH+'data/'+x for x in cell_types], sys.argv[1])
    
    
