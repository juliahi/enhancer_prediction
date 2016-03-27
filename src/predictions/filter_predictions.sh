#best N predictions
N=10000
  
for tissue in 'brain' 'heart'
do
    python -m src.predictions.filter_predictions "${DB}_${tissue}__4mers" "predicted_${DB}_${tissue}__4mers_$N/all.bed"  best $N
    python -m src.predictions.filter_predictions "${DB}_${tissue}_h1hesc_4mers" "predicted_${DB}_${tissue}_h1hesc_4mers_$N/all.bed"  best $N
done

#predictions over 80%
TH=${CUTOFF}
for tissue in 'brain' 'heart'
do
    python -m src.predictions.filter_predictions "${DB}_${tissue}__4mers" "predicted_${DB}_${tissue}__4mers_$TH/all.bed" threshold $TH
    python -m src.predictions.filter_predictions "${DB}_${tissue}_h1hesc.list_4mers" "predicted_${DB}_${tissue}_h1hesc.list_4mers_$TH/all.bed"  threshold $TH
done
