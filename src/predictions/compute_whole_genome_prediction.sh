
for chrom in {1..22} #X
do
    echo predicting chromosome $chrom 
    for tissue in 'heart' 'brain' 'both'
    do
         python -m src.predictions.classify "predictions/$WINDOW/chr$chrom.id" $DB $tissue "-" 4mers &
#          python -m src.predictions.classify "predictions/$WINDOW/chr$chrom.id" $DB $tissue "h1hesc" 4mers &
#          python -m src.predictions.classify "predictions/$WINDOW/chr$chrom.id" $DB $tissue "tier12" 4mers &
    done
done
wait

