import os

DATAPATH=os.environ["DATAPATH"]


RESULTSPATH=os.environ["RESULTSPATH"]

IMGPATH=RESULTSPATH+"images/"
SUMMARYFILE=RESULTSPATH+"summary.csv"

USEGC = False

#http://www.ncbi.nlm.nih.gov/projects/genome/assemly/grc/human/data/index.shtml
#GRCh37
chr_lens = dict([('chr1', 2000621), #('chr1', 249250621), 
		('chr2',  243199373), ('chr3',  198022430),
                 ('chr4',  191154276), ('chr5',  180915260), ('chr6',  171115067), 
                 ('chr7',  159138663), ('chr8',  146364022),('chr9',  141213431), 
                 ('chr10', 135534747), ('chr11', 135006516), ('chr12', 133851895),
                 ('chr13', 115169878), ('chr14', 107349540), ('chr15', 102531392), 
                 ('chr16', 90354753), ('chr17', 81195210), ('chr18', 78077248), 
                 ('chr19', 59128983), ('chr20', 63025520), ('chr21', 48129895), 
                 ('chr22', 51304566), ('chrX',  155270560), ('chrY',  59373566)]) 

data_lists = ["h1hesc.list"
        #, "tier1.list"
        , "tier12.list"
        ]

cell_types = ['Cd20bis', 'Cd20ro01778bis',  'Cd20ro01794bis',  'Gm12878',  
    'H1hesc',  
    'Huvec',  'Huvecbis',  'K562',  'Monocd14ro1746bis'
    ]

cv_folds = 10 #x-fold crossvalidation
N_trees = 100 #random forest parameter
N_repeats = 5 #repeats to count auc

N_pools = 6 #n of threads

GENOME_FILE = "female.hg19.fa"


TSS_FILE = DATAPATH+'ensemble_hg19_tss_sorted.bed'

#whole genome predictions parameters
WINDOW=int(os.environ["WINDOW"])
THRESHOLD=os.environ["CUTOFF"]
STEP=WINDOW/2
WIGGLENAME="predictions/%d/"%WINDOW
WIGGLEPATH=DATAPATH+WIGGLENAME
VISTA_FILE=os.environ["DB"]

def get_name(filename):
    k = filename[:filename.find(".")]
    k = k.split('/')[-1].replace('wgEncode', '').replace('Histone', '').replace('StdAln', '').replace('Reps', '').replace('RawRep', '_').replace('Sig', '').replace('Std', '')
    k = k.replace('H3k', '_H3k').replace('H4k', '_H4k').replace('H2az', '_H2az') 
    return k




