#!/bin/bash
## Genome Size Pipeline ##
## Given a histogram file from flowJo and post Excel, returns the expected estimated genome size for D. mel ##

# DIRECTORIES
mainDir="/Users/alisonhanh/Desktop/data/genome_size/pipeline/"
outDir=${mainDir}output/
pdfDir=${outDir}pdf/

# EXCEL  output - number after "line" represents day batches
# Excel's Macros VBA script will place the output files in this directory with prefixes indicating the day batch
xlsDir="/Volumes/ALISONN/facs/Y-lines2016/flowJo_files/"
fixedTsvDir=${mainDir}tsvFiles/
histTsvDir=${mainDir}histTsv/

# Genome Estimate files to handle specified in this txt
fileList=${mainDir}pipeline_input.txt

## PHASE 1: FIX TSV FROM ^M TO NEWLINE IN UNIX ##
## PHASE 1 STATUS: PASS ##
fixTsvFormat () {
	# the directories are fixed in the secondary bashscript so you have to modify it and then call
	bashScript=${mainDir}fix_tsv_format.sh
	$bashScript
} 

## PHASE 2: CONVERSION FROM TSV TO R HISTOGRAM COMPATIBLE FORMAT IN PYTHON ##
## PHASE 2 STATUS: PASS ##
convertTsv() {
	pythonScript=${mainDir}tsv_to_hist_R.py
	## for loop through all of the files in this folder - you can modify which batch you loop over
	cd $fixedTsvDir
	for file in *.tsv; do
		$pythonScript $file ${histTsvDir}$file
		mv $file finished
	done
	cd $mainDir
}

## PHASE 3: GENERATION OF HISTOGRAM AND PEAK ESTIMATION IN R ##
getRPlot() {
	# in the histTsvDir directory, loop through all of the tsv and make the plots
	plotScript=${mainDir}calculate_genome_size.R
	cd $histTsvDir
	for file in *.tsv; do
		prefix=$(echo $file | cut -f 1 -d '.')
		outTsv=${outDir}yLine_estimates_master.tsv
		outPdf=${pdfDir}$prefix.pdf

		# this script will append the results of every estimate into the master sheet so you can read in the .tsv for easy anova testing
		# first argument is the hist_input.tsv, second is the master_output.tsv, third is the pdf file
		Rscript $plotScript $file $outTsv $outPdf
		mv $file finished
	done
	cd $mainDir
}

## PHASE 4: CALL ALL FUNCTIONS ##

# fixTsvFormat
# convertTsv
# getRPlot

## ALL PHASES DONE ##
