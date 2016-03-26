
N=$WINDOW

for name in $DB randoms_heart randoms_brain 
do
	python -m src.promoters.vista_window $N $name$N $name

	python -m src.preprocessing.histmods $name$N &
	python -m src.preprocessing.gen_kmers $name$N 4 &
	wait
done

