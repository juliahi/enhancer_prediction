import os

from ..shared import *
from copy import copy
import promoters

class Region:
    def __init__(self, chromosome, start, end, value=0, id=''):
        self.chromosome = chromosome
        self.start = int(start) #1-based
        self.end = int(end)     #including (to match Position class with enhancers)
        self.value = float(value)
        self.seq = ""
        self.id = id
    def __repr__(self):
       return "%s\t%d\t%d\t.\t%f\t." % (self.chromosome, self.start-1, self.end, self.value)  #3.12.15 added -1 -> 0-based, excluding
    def to_fasta(self):
        s = ">%s:%s-%s | %s" % (self.chromosome, self.start, self.end, self.id) #1-based
        if self.seq != '':
            s += '\n' + str(self.seq)
        return s
    

def overlaps(reg1, reg2):
    if reg2[0] < reg1[0] < reg2[1]:
        return True
    if reg2[0] < reg1[1] < reg2[1]:
        return True
    if reg1[0] < reg2[0] < reg1[1]:
        return True
    return False
    
def get_scores(regions, wiggle):
    #regions - list of [chrom, start, end]
    #wiggle - dict[chromosome: [Regions]]
    scores = []
    for chrom, start, end in regions:
        values = []
        if chrom in wiggle:
            for region in wiggle[chrom]:
                if overlaps((region.start, region.end), (start, end)):
                    values.append( region.value )
                else:
                    if values != []:
                        break
            if len(values) > 0:
                scores.append(sum(values)/len(values))
    return scores
        



def read_wiggle(filename):
    #Returns dictionary with chromosomes and lists of Regions
    with open(filename) as f:
        f.readline()
        
        result = {}
        while True:
            line = f.readline()
            if not line: break
            if line.startswith('variableStep'):
                span = int(line.split()[2].split('=')[1])
                chrom = line.split()[1].split('=')[1].split('.')[0]
                result[chrom] = []
            elif len(line.split()) == 2:
                result[chrom].append(Region(chrom, int(line.split()[0]), int(line.split()[0])+span-1, line.split()[1]))
                
        return result

def read_all(name):
    #read all wig files from directory name
    wiggle = {}
    for i in range(1, 23)+['X']:
        MNAME = WIGGLEPATH + name + ('/chr%s.id.wig'%str(i))
        print MNAME, os.path.exists(MNAME)
        if os.path.exists(MNAME):        
            wiggle.update(read_wiggle(MNAME))
    return wiggle

def read_some(name, chrs):
    wiggle = {}
    for i in chrs: 
        MNAME = WIGGLEPATH + name + ('/chr%s.id.wig'%str(i))
        wiggle.update(read_wiggle(MNAME))
    return wiggle

def read_bed(name):
    result=[]
    with open(name, 'r') as f:
        for line in f:
            
            prom = line.strip().split('\t')
            if len(prom) < 4 or line[0] == '#':
                continue
            #if len(prom) < 6:
                #strand = '.'
            #elif prom[5] == '+' or prom[5] == "1":
                #strand = '+'
            #else:
                #strand = '-'
            
            if prom[0] == 'MT': prom[0] = 'M'
            if prom[0][:3] != 'chr':
                prom[0] = 'chr'+prom[0]
            #print prom[0], prom[1], prom[2], prom[3], prom[4]
            if len(prom) < 6: #not a proper bed file (like Erwin)
                name = prom[1]
                score = float(prom[3].split(';')[0])
            else:
                name = prom[3]
                score = float(prom[4])
            result.append(Region(prom[0], int(prom[1])+1, int(prom[2]), score, name))
    return result



def read_bed_dict(name):
    result={}
    with open(name, 'r') as f:
        for line in f:
            
            prom = line.strip().split('\t')
            if len(prom) < 4 or line[0] == '#':
                continue
            #if len(prom) < 6:
                #strand = '.'
            #elif prom[5] == '+' or prom[5] == "1":
                #strand = '+'
            #else:
                #strand = '-'
            
            if prom[0] == 'MT': prom[0] = 'M'
            if prom[0][:3] != 'chr':
                prom[0] = 'chr'+prom[0]
            #print prom[0], prom[1], prom[2], prom[3], prom[4]
            if len(prom) < 6: #not a proper bed file (like Erwin)
                name = prom[1]
                score = float(prom[3].split(';')[0])
            else:
                name = prom[3]
                score = float(prom[4])
            if prom[0] in result:
                result[prom[0]].append(Region(prom[0], int(prom[1])+1, int(prom[2]), score, name))
            else:
                result[prom[0]] = [Region(prom[0], int(prom[1])+1, int(prom[2]), score, name)]
                
    return result