


db="${DB}1500 promoters"
for tissue in 'both' 
do 
    positives="randoms_both"
    negatives="predicted_promoters_4mers_balanced_0.8"
    python -m src.rf.run --db $db -p $positives -n $negatives                     --kmers 4mers       --distinct    -o "classifiers/class_balancedpromoters_${positives}_4mers"    --saveclass     &
    
done


wait
