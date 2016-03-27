

python -m src.predictions.cut_seq $WINDOW

for chrom in $CHROMOSOMES
do
    echo $chrom
    python -m src.preprocessing.gen_kmers "predictions/$WINDOW/chr$chrom.id" 4 &
done

wait

for chrom in $CHROMOSOMES
do
    echo $chrom
    python -m src.preprocessing.histmods "predictions/$WINDOW/chr$chrom.id"  &
    sleep 100
done

wait


