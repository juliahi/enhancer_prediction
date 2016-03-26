saveclass='--saveclass'

for i in {1..5} 
do
    negative='predicted_promoters'
    for tissue in 'heart' 'brain' 'both' 
    do 
        args=" --db $db promoters -p ${tissue}_notss --distinct $saveclass"
        python -m src.rf.run $args -n ${negative}_4mers_${tissue}_0.8          --kmers 4mers       -o "promoters_${DB}_${tissue}_4mers"  &
        python -m src.rf.run $args -n ${negative}_tier12_4mers_${tissue}_0.8 --histmods tier12.list               -o "promoters_${DB}_${tissue}_tier12" &
        python -m src.rf.run $args -n${negative}_tier12_4mers_${tissue}_0.8 --histmods tier12.list  --kmers 4mers -o "promoters_${DB}_${tissue}_tier12_4mers"  &
    
    done
done
wait


