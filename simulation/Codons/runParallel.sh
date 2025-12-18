#!/bin/bash

nthreads=4

rm templates/rep*/*.log

for ((i=1;i<=nthreads;i++)); do

	echo "$i/$nthreads"
	nohup bash runParallel_b.sh $i $nthreads & disown
	#bash runParallel_b.sh $i $nthreads

done

