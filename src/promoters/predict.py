    
import os
from Bio import SeqIO

from ..shared import *
from copy import copy
import promoters
from read_wiggle import *
    

    
def save_promoters_predictions(tss_dict, classifier, cutoff, outname):
    wiggle = read_wiggle(WIGGLEPATH+classifier)
    print 'Finding promoter predictions', classifier

    if '/' in outname: 
        dirname = os.path.dirname(outname)
        if not os.path.exists(DATAPATH+dirname):
            print 'creating directory', DATAPATH+dirname
            try:
                os.makedirs(DATAPATH+dirname)
            except:
                pass
    outfile = open(DATAPATH+outname+".tmp", 'w+')
    count = 0
    countn = 0
    
    
    for chrom, regions in wiggle.iteritems():
        tss_list = tss_dict[chrom]
        for region in regions:
            for start, strand in tss_list:
                if strand == 1 or strand == '+':
                    window_start = start - start%STEP - STEP
                else:
                    window_start = start - start%STEP
                if region.start-1 == window_start:
                        if region.value > cutoff:
                            name = (region.start-1) / STEP   #id is 
                            outfile.write('%s\t%d\t%d\telement_%d_tss_%d\t%f\t%s\n' % (chrom, region.start-1, region.end, name, start, region.value, strand))
                            count += 1
                        countn += 1
                        break
        outfile.flush()
        print 'chromosome %s finished' % chrom
    print "%d promoters from %d over cutoff %f" % (count, countn, cutoff)
    outfile.close()
    os.system("cat %s.tmp | sort -k 1,1 -k 2,3n -u > %s_sorted.bed"%(DATAPATH+outname, DATAPATH+outname))
    os.system("rm %s.tmp"%(DATAPATH+outname))
    


    
def bed_to_fasta_seq(name, outname):
    enhancers = promoters.load_dict(name)
    f = open(DATAPATH+outname, 'w')
    
    for record in SeqIO.parse(DATAPATH+GENOME_FILE, "fasta"):
        chromseq = record.seq
        if record.id in enhancers:
            for e in enhancers[record.id]:
                seq = str(chromseq[e.start:(e.end+1)])
                f.write(e.to_fasta()+'\n'+seq+"\n")
            
        print 'chromosome %s finished'%record.id
    f.close()


if __name__ == '__main__':
    import sys
    print sys.argv
    if len(sys.argv) <= 3:
        print "USAGE: predict.py input_filename output_name promoters|sequences [threshold]"
        sys.exit(1)
    f = DATAPATH + sys.argv[1]
    
    
    if len(sys.argv) < 5:
        threshold = 0
    else:
        threshold = float(sys.argv[4])
    if sys.argv[3] == 'promoters':
        tss_list = promoters.read_tss()
        save_promoters_predictions(tss_list, sys.argv[1], threshold, sys.argv[2])
    if sys.argv[3] == 'sequences':
        print "bed to fasta"
        bed_to_fasta_seq(sys.argv[1], sys.argv[2])
        
    