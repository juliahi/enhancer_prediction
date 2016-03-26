from Bio import SeqIO
import re
from ..shared import *

class Position:
    def __init__(self, chromosome, start, end, id, positive, tissue, seq=''):
        self.chromosome = chromosome
        self.start = int(start)   #1-based
        self.end = int(end)  #including!!!!!!!!! TODO
        self.id = int(id)
        self.positive = positive
        self.tissue = tissue
        self.seq = seq
    def __str__(self):
        if self.positive:
            strtiss = " | ".join(self.tissue)
            s = ">Human| %s:%s-%s | element %s | positive | %s" % (self.chromosome, self.start, self.end, self.id, strtiss)
        else:
            s = ">Human| %s:%s-%s | element %s | negative" % (self.chromosome, self.start, self.end, self.id)
        if self.seq != '':
            s += '\n' + str(self.seq)
        return s 
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
    def __hash__(self):
        """Override the default hash behavior (that returns the id or the object)"""
        return hash(tuple(sorted(self.__dict__.items())))

    
    
    def name(self):
        if self.positive:
            return "%d_positive_%s"% (self.id, self.tissue)
        else:
            return "%d_negative_%s" % (self.id, self.tissue)
    def is_neural(self):
        neural_tissues = ["neural tube", "cranial nerve", "hindbrain (rhombencephalon)", "forebrain", "midbrain (mesencephalon)", "brain", "neural"]
        found = [s for s in neural_tissues if s in self.tissue]
        if found == []:
            return False
        else:
            return True
    def is_limb(self):
        if "limb" in self.tissue:
            return True
        else:
            return False
    def is_heart(self):
        if "heart" in self.tissue:
            return True
        else:
            return False
    def is_random(self):
        #return (int(self.id) >= 10000)
        return ("random" in self.tissue)

tissue_types = ["heart", "neural", "limb", "positives", "vista"]



def load_enhancers(filename, tissues=[], withseq = False):
    handle = open(DATAPATH+filename, "rU")
    positions = []
    for record in SeqIO.parse(handle, "fasta") :
        #print record
        rec = record.description.split('|')
        tissue = tissues[:]
        if len(rec) < 3: #load FANTOM
            pos = rec[-1].strip().split(':')
            startend = pos[1].split('-')
            positive = True
            if len(pos) > 2:
                id = pos[2]
            else:
                id = -1
        else:
            pos = rec[1].strip().split(':')
            startend = pos[1].split('-')
            id = rec[2].strip()[len("element "):]
            if rec[3].strip().startswith("positive"):
                positive = True
            elif rec[3].strip().startswith("negative"):
                positive = False
            else: #not defined - predictions?
                positive = False
            for i in xrange(4, len(rec)):
                tissue.append(re.sub('\[[0-9\/]*\]', '', rec[i]).strip())
        if withseq:
            positions.append(Position(pos[0], startend[0], startend[1], id, positive, tissue, record.seq))
        else:
            positions.append(Position(pos[0], startend[0], startend[1], id, positive, tissue))
    handle.close()
    return positions



def load_enhancers_with_seq(filename, tissues=[]):
    return load_enhancers(filename, tissues, True)



def filter_type(enhancers, typ):
    if typ == "heart":
        return [e for e in enhancers if e.is_heart()]
    if typ == "neural":
        return [e for e in enhancers if e.is_neural()]
    if typ == "limb":
        return [e for e in enhancers if e.is_limb()]
    if typ == "positives":
        return [e for e in enhancers if e.positive]
    if typ == "vista":
        return [e for e in enhancers if not e.is_random()]
    return enhancers


def save_enhancers(enhancers, filename):
    f = open(DATAPATH+filename, 'w')
    for e in enhancers:
       f.write(str(e)+"\n")
   



def load_bed(filename, positive=False, tissues=[], startid=1):
    handle = open(DATAPATH+filename, "r")
    positions = []
    id = startid
    for line in handle:
        info = line.split('\t')
        positions.append(Position(info[0], int(info[1])+1, int(info[2]), id, positive, tissues))
        id += 1
    return positions

def save_to_bed(enhancers, filename):
    f = open(DATAPATH+filename, 'w')
    for e in enhancers:
        name = "%d_%s_%s" % (e.id, e.positive, e.tissue)
        f.write(e.chromosome + "\t" + str(e.start-1) + "\t" + (str(e.end)) + "\t" + name + "0\t.\n")
        # changed 24.09.15 
    f.close() 


def add_seqs(enhancers, output_file):                            #add sequences
    record_dict = SeqIO.index(DATAPATH+GENOME_FILE, "fasta")
    for enh in enhancers:
        if enh.seq == '':
            enh.seq = str(record_dict[enh.chromosome].seq[(enh.start-1):(enh.end)])  
            #cause position is 1-based, including end
    save_enhancers(enhancers, output_file) 




