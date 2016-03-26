

WINDOW=1500
# for chrom in {1..22} #X
# do
#     echo predict $chrom 
#     for tissue in  'heart' 'brain'
#     do
#          python -m src.predictions.classify_two_step "predictions/$WINDOW/chr$chrom.id" vista $tissue "-" 4mers &
# #          python -m src.predictions.classify "predictions/$WINDOW/chr$chrom.id" vista $tissue "h1hesc" 4mers &
# #          python -m src.predictions.classify "predictions/$WINDOW/chr$chrom.id" vista $tissue "tier12" 4mers &
#     done
# done
# wait
# 


### only second step
for chrom in 1 #{1..22} #X
do
    echo predict $chrom 
#     python -m src.predictions.classify "predictions/$WINDOW/chr$chrom.id" 'balancedpromoters' 'randoms_both' "-" 4mers &
    python -m src.predictions.classify "predictions/$WINDOW/chr$chrom.id" 'vista' 'heart_notss' "-" 4mers &
done
wait