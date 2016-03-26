boruta='--boruta'
for i in {1..5}
do
    db="$DB1500 promoters"
    tissue='both_notss'
    tissue_neg='predicted_promoters_4mers_both_0.8'
    
    python -m src.rf.run --db $db -p $tissue -n $tissue_neg                    --kmers 4mers       --distinct    -o "both_enhancers_vs_promoters_4mers" $boruta &
    python -m src.rf.run --db $db -p $tissue -n $tissue_neg --histmods tier12.list                  --distinct    -o "both_enhancers_vs_promoters_tier12" $boruta &
    python -m src.rf.run --db $db -p $tissue -n $tissue_neg --histmods tier12.list  --kmers 4mers    --distinct   -o "both_enhancers_vs_promoters_tier12_4mers" $boruta &
    boruta=''
    sleep 10
done
wait


