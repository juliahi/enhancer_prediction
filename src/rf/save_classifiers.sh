db=$DB
for tissue in 'heart' 'brain' 'both'; do 
        python -m src.rf.run --db $db -p $tissue -n randoms_$tissue                     --kmers 4mers       --distinct    -o "classifiers/class_${db}_${tissue}_4mers"    --saveclass     &
        python -m src.rf.run --db $db -p $tissue -n randoms_$tissue --histmods h1hesc.list                  --distinct    -o "classifiers/class_${db}_${tissue}_h1hesc"   --saveclass     &
        python -m src.rf.run --db $db -p $tissue -n randoms_$tissue --histmods h1hesc.list  --kmers 4mers    --distinct   -o "classifiers/class_${db}_${tissue}_h1hesc_4mers" --saveclass &
        python -m src.rf.run --db $db -p $tissue -n randoms_$tissue --histmods tier12.list                  --distinct    -o "classifiers/class_${db}_${tissue}_tier12"   --saveclass     &
        python -m src.rf.run --db $db -p $tissue -n randoms_$tissue --histmods tier12.list  --kmers 4mers   --distinct    -o "classifiers/class_${db}_${tissue}_tier12_4mers" --saveclass &
done
wait
