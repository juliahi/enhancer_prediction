boruta='--boruta'
for i in {1..$REPEATS }; do
    for tissue in 'heart' 'brain' 'both'; do 
            python -m src.rf.run --db $DB -p $tissue -n randoms_$tissue                     --kmers 4mers       --distinct    -o "${DB}_${tissue}_4mers" $boruta &
            python -m src.rf.run --db $DB -p $tissue -n randoms_$tissue --histmods h1hesc.list                  --distinct    -o "${DB}_${tissue}_h1hesc" $boruta &
            python -m src.rf.run --db $DB -p $tissue -n randoms_$tissue --histmods h1hesc.list  --kmers 4mers    --distinct   -o "${DB}_${tissue}_h1hesc_4mers" $boruta &
            python -m src.rf.run --db $DB -p $tissue -n randoms_$tissue --histmods tier12.list                  --distinct    -o "${DB}_${tissue}_tier12" $boruta &
            python -m src.rf.run --db $DB -p $tissue -n randoms_$tissue --histmods tier12.list  --kmers 4mers   --distinct    -o "${DB}_${tissue}_tier12_4mers" $boruta &
	    wait
    done
    boruta=''
done
