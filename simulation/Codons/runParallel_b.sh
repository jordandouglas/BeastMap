#!/bin/bash

END=200


cd templates



for ((i=1;i<=END;i++)); do

	j=$(($i-$1))
	#echo $j
	if [ `expr $j % $2` -eq 0 ]
	then
	    echo $i
		cd rep$i
		~/beast/bin/beast -overwrite -df var.seq.json ../../run.xml
		cd ../
	fi


done


cd ../


