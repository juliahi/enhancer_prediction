echo Usage ./run_all.sh DATAPATH RESULTSPATH

if [ $# -ne 2 ]; then
	exit 1
fi

DATAPATH=$1/
DB=vista
RESULTSPATH=$2/
export DATAPATH DB RESULTSPATH
WINDOW=1500
export WINDOW
REPEATS=1
export REPEATS
CUTOFF=0.1
export CUTOFF

echo "Preprocessing"
#python -m src.preprocessing.histmods $DB
#python -m src.preprocessing.gen_kmers $DB 4

#random sequences
#python -m src.preprocessing.random_seq gen $DB heart randoms_heart 100000    #100000 = starting sequence id
#python -m src.preprocessing.random_seq gen $DB brain randoms_brain 200000

#for tissue in 'heart' 'brain'
#do
#    python -m src.preprocessing.histmods randoms_$tissue 
#    python -m src.preprocessing.gen_kmers randoms_$tissue 4 
#done



#########################
#Training:
#Cross-validation:
mkdir -p $RESULTSPATH
mkdir -p $RESULTSPATH/datasets
mkdir -p $RESULTSPATH/classifiers
echo "DB	POS	NEG	kmers	hmods	dist	ntrees	use_GC	cv_folds	used_class	AUC	date	cut_at	output_dir	POS_size	NEG_size" > $RESULTSPATH/summary.csv

echo "Training"
#./src/rf/run_crossval.sh

#train and save
#./src/rf/save_classifiers.sh

#compute pvalues of changes of AUC
#python -m src.statistics.pvalues



############################
#whole genome predictions
echo "Whole genome predictions"
mkdir -p "${RESULTSPATH}predictions/${WINDOW}"
CHROMOSOMES="1" #"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"
export CHROMOSOMES
# generate sequences with window 1500 overlap 750, without N in sequence, count features 
#./src/predictions/prepare_data.sh 

#predict, write in wig format to RESULTSPATH/predictions/chrN.id.wig
#./src/predictions/compute_whole_genome_prediction.sh

#Select best predictions:
#./src/predictions/filter_predictions.sh

############################
#Promoters 
echo 'Promoters recognition'

#Prepare training set -- adjust vista and randoms lengths to 1500
#./src/promoters/vista_window.sh

#./src/promoters/select_active_promoters.sh 

#classify
#./src/promoters/run_crossval.sh
./src/promoters/save_class.sh


############################
echo 'Two-step classifier'
#Train second step classifiers
./src/rf/save_classifiers_second_step.sh

#Whole genome predictions
./src/predictions/compute_whole_genome_prediction_two_step.sh







###########################
#Plots 

#./src/paper_plots.sh



