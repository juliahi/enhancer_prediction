saveclass='--saveclass'

for i in {1..$REPEATS } 
do
    negative='predicted_promoters'
    for tissue in 'heart' 'brain' 'both' 
    do 
        args=" --db $DB promoters -p ${tissue}_notss --distinct $saveclass"
        python -m src.rf.run $args -n ${negative}_4mers_${tissue}_${CUTOFF}          --kmers 4mers       -o "promoters_${DB}_${tissue}_4mers"  &
        #python -m src.rf.run $args -n ${negative}_h1hesc_4mers_${tissue}_${CUTOFF} --histmods h1hesc.list               -o "promoters_${DB}_${tissue}_h1hesc" &
        #python -m src.rf.run $args -n${negative}_h1hesc_4mers_${tissue}_${CUTOFF} --histmods h1hesc.list  --kmers 4mers -o "promoters_${DB}_${tissue}_h1hesc_4mers"  &
    
    done
done
wait


