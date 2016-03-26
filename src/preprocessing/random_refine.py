from load_vista import load_enhancers, Position, save_enhancers
import numpy
#import matplotlib.pyplot as plt
import os
from scipy import stats
import random
from Bio import SeqIO
import histonemod
import histstats
from lens import chr_lens

home = os.getenv("HOME") + "/"

def change_random(randoms_file, to_remove):
    enhancers = load_enhancers(randoms_file)
    record_dict = SeqIO.index(home+"data/female.hg19.fa", "fasta")
    to_return = []
    for e in enhancers:
        if e.id not in to_remove:
            continue
        while True: 
            l = e.end - e.start
            start = random.randint(0,chr_lens[e.chromosome]-l-1)
            if record_dict[e.chromosome].seq[start:(start+l)].find('N') == -1:
                to_return.append(Position(e.chromosome, start, start + l, e.id, False, []))
                break
            else:
                print e.chromosome
                
    save_enhancers(to_return, randoms_file+".change")
    return to_return
    
    
def replace_randoms(randoms_old, randoms_replace, randoms_new):
    rold = load_enhancers(randoms_old)
    rreplace = load_enhancers(randoms_replace)
    rnew = []
    if rreplace == []:
        return
    e2 = rreplace.pop(0)
    for e in rold:
       if e.id == e2.id:
           rnew.append(e2)
           if rreplace != []:
               e2 = rreplace.pop(0)
       else:
           rnew.append(e)
    
    save_enhancers(rnew, randoms_new)
   

def replace_results(randoms_old, randoms_replace, directory, suffix, replace_suffix, new_suffix):
    rold = load_enhancers(randoms_old)
    rreplace = load_enhancers(randoms_replace)
    if rreplace == []:
        return
    e2 = rreplace.pop(0)
    lines = []
    for e in rold:
       if e.id == e2.id:
           lines.append(True)
           if rreplace != []:
               e2 = rreplace.pop(0)
       else:
           lines.append(False)
    for filename in glob(directory + "/*" + suffix):
       frep = open(filename + replace_suffix)
       fnew = open(filename + new_suffix)
       for i, line in enumerate(open(filename)):
           if lines[i]:
              fnew.write(frep.readline()) 
           else:
              fnew.write(line)
       frep.close()
       fnew.close()
       


def refine(name, plot_img):
    filename = home + "data/" + name
    to_change = histstats.stats(filename, name, plot_img)
    if to_change == []:
        return
    print "changing %d sequences" % len(to_change)
    change_random(filename, to_change)
    replace_randoms(filename, filename+".change", filename)
    histonemod.change_all(filename, filename+".change", name)
    refine(name, False)



home = os.getenv("HOME") + "/"
if  __name__ =='__main__':
    #refine("randoms_heart", True)
    #refine("randoms_limb", True)
    #refine("randoms_neural", True)
    refine("randoms_positives", True)
    #refine("randoms_vista", True)




