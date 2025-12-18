#!/bin/bash

END=100

cd templates



for ((i=1;i<=END;i++)); do

	echo $i
	cd rep$i
	~/beast/bin/beast -overwrite -df var.seq.json ../../run.xml
	cd ../


done





