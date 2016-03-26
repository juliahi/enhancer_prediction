

for tissue in 'heart' 'brain' 'both' 
do
    for chr in {1..22} #X
    do
         python -m src.promoters.predict "vista_${tissue}__4mers/chr${chr}.id.wig" "predicted_promoters_4mers_${tissue}_0.8/chr${chr}" promoters 0.8 &
          python -m src.promoters.predict "vista_${tissue}_h1hesc.list_4mers/chr${chr}.id.wig" "predicted_promoters_h1hesc_4mers_${tissue}_0.8/chr${chr}" promoters 0.8 &
          python -m src.promoters.predict "vista_${tissue}_tier12.list_4mers/chr${chr}.id.wig" "predicted_promoters_tier12_4mers_${tissue}_0.8/chr${chr}" promoters 0.8 &
    done
    wait
done

