#!/bin/bash

# @param: excel sheet --> tab-delimited file for processing genome size estimates
# @return: a tab-delimited file without ^M as newline characters for unix/bash compatibility
# instead of using sed for whatever reason the special char ^M won't work
# use tr function instead 

currDir="/Volumes/ALISONN/facs/Y-lines2016/flowJo_files/"
outDir="/Users/alisonhanh/Desktop/data/genome_size/pipeline/tsvFiles/"
cd $currDir

for file in *.tsv; do
	# this gets the STRAIN.tsv portion of the excel-to-tsv file
	sample=${outDir}$file
	echo "This is the output and its location"
	echo $sample
	echo \n
	tr '\r' '\n' < $file > $sample
	# eventually I'll kill the old .tsv files with the ^M special char
	# rm $file
done
