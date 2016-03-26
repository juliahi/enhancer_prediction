echo Usage ./run_all.sh DATAPATH name_of_db RESULTSPATH

if [ $# -ne 3 ]; then
	exit 1
fi

DATAPATH=$1/
DB=$2
RESULTSPATH=$3/
export DATAPATH DB RESULTSPATH
WINDOW=1500
export WINDOW

python -m src.preprocessing.histmods $DB
python -m src.preprocessing.gen_kmers $DB 4

#random sequences
python -m src.preprocessing.random_seq gen $DB heart randoms_heart 100000    #100000 = starting sequence id
python -m src.preprocessing.random_seq gen $DB brain randoms_brain 200000

for tissue in 'heart' 'brain'
do
    python -m src.preprocessing.histmods randoms_$tissue 
    python -m src.preprocessing.gen_kmers randoms_$tissue 4 
done



#########################
#Training:
#Cross-validation:
mkdir $RESULTSPATH
touch $RESULTSPATH/summary.csv

./src/rf/run_crossval.sh


#train and save
./src/rf/save_classifiers.sh

#compute pvalues of changes of AUC
python -m src.statistics.pvalues


############################
#whole genome predictions
# generate sequences with window 1500 overlap 750, without N in sequence, count features 
./src/predictions/prepare_data.sh 

#predict, write in wig format to RESULTSPATH/predictions/chrN.id.wig
./src/predictions/compute_whole_genome_predictions.sh

#Select best 1000 predictions:
./src/predictions/filter_predictions.sh

############################
#Promoters 


#Prepare training set -- adjust vista and randoms lengths to 1500
./src/promoters/vista_window.sh

./src/promoters/select_active_promoters.sh 

#classify
./src/promoters/run_crossval.sh
./src/promoters/save_class.sh

############################

#Train second step classifiers
./src/rf/save_classifiers_second_step.sh

#Whole genome predictions
./src/predictions/compute_whole_genome_prediction_two_step.sh







###########################
#Plots 

./src/paper_plots.sh



