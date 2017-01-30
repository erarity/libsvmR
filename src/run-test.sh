#!/bin/bash

vlens="500" # 2000 5000"
sparses="0.50 0.75 0.90"
skips="0.10 0.25 0.50"
rows="500"
reps="1"

for vlen in $vlens; do
    for sparse in $sparses; do
	for skip in $skips; do
	    for row in $rows; do
		for rep in $reps; do
		    # echo Running with params $vlen $sparse $skip $row $rep
		    ./a.out $vlen $sparse $skip $row $rep 
		done
	    done
	done
    done
done | tee all-output.txt


		    
	    
