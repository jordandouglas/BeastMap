#!/bin/bash



# Prepare templates to do WCSS
rm -r templates
mkdir -p templates

# Sample from prior
mkdir -p truth
cd truth
~/beast/bin/beast -overwrite ../truth.xml
cd ../



# Put sampled variables into json files
Rscript scripts/populate.R



cd templates
for rep in rep*/;
do
	cd $rep
	~/beast/bin/beast -overwrite -df var.json ../../seqsim.xml


	echo "---------------"
	echo $rep
	echo "---------------"

	# Put sequences into json file
	Rscript ../../scripts/putSequencesIntoJson.R


	cd ../

done
cd ../

bash runParallel.sh
# #bash run.sh







