#!/bin/bash
# Job name:
#SBATCH --job-name=polish
#SBATCH --account=co_bachtrog
#SBATCH --partition=savio2_bigmem
#SBATCH --qos=bachtrog_bigmem2_normal
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=120:00:00

### Polish genome assemblies from Canu assembler 3x with RACON and 1x with PILON

## Command(s) to run:
export PATH=$PATH:/global/home/groups/co_bachtrog/programs/minimap2/minimap2-master/
export PATH=$PATH:/global/home/groups/co_bachtrog/programs/miniasm
export PATH=$PATH:/global/home/groups/co_bachtrog/programs/samtools-1.5
export PATH=$PATH:/global/home/groups/co_bachtrog/programs/bwa_version_0.7.15

runRacon () {

	## Command(s) to run:
	export PATH=$PATH:/global/home/groups/co_bachtrog/programs/minimap2/minimap2-master/
	export PATH=$PATH:/global/home/groups/co_bachtrog/programs/miniasm
	export PATH=$PATH:/global/home/groups/co_bachtrog/programs/samtools-1.5
	#export PATH=$PATH:/global/home/groups/co_bachtrog/programs/racon/bin
	export PATH=$PATH:/global/scratch/dmai/software/racon/build/bin
	export PATH=$PATH:/global/home/groups/co_bachtrog/programs/bwa_version_0.7.15
	##export PATH=$PATH:/global/home/users/bracewel/Programs/samtools-0.1.19

	mapDir="/global/scratch/alisonn/02_mapped"

	refDir=$1
	readsDir=$3

	genome=${refDir}/$2
	genomePrefix=$( echo $2 | cut -f 1 -d '.' )
	longReads=${readsDir}/$4
	#longReads=$3
	iteration=$5

	output=${refDir}/${genomePrefix}.reads_mapped.racon_${iteration}.paf
	updatedGenome=${refDir}/${genomePrefix}.racon_${iteration}.contigs.fasta

	echo "Running RACON to polish with long reads"
	echo "genome:        $genome"
	echo "reads:         $longReads"
	echo "iteration:     $iteration"
	echo "mapped output: $output"
	echo "new genome:    $updatedGenome"

	ls $genome
	ls $longReads
	ls $output	

	echo "starting minimap2"
	minimap2 $genome $longReads > $output
	echo "starting racon"
	#racon
	racon -t 24 -u $longReads $output $genome > $updatedGenome 
	## if you get the error chunk size too small, make sure your fastq.gz file has quality scores.. this is usually the issue

	### old racon version command is below
	###racon -t 24 $longReads $output $genome $updatedGenome 

	echo " "; echo " "

}

export -f runRacon

runPilon () {

	genomeDir=$1
	readsDir=$3
	mappedDir="/global/scratch/alisonn/02_mapped/dpse-nanopore"

	# genome reference is in this format: dpse124Y_nanopore_canu_v1.racon_3.contigs.fasta
	prefix=$( echo $2 | sed -E 's/(\w+.racon_[0-9]{1,12}).contigs.fasta/\1/g' )
	genome=${genomeDir}/$2
	read1=${readsDir}/$4
	read2=${readsDir}/$5
	iteration=$6
	samOut=${mappedDir}/${prefix}.illuminaMap.sam
	bamOut=${mappedDir}/${prefix}.illuminaMap.sorted.bam
	updatedGenome=${genomeDir}/${prefix}.pilon_${iteration}.contigs.fasta

	echo "Running PILON to polish with short reads"	
	echo "prefix:          			$prefix"
	echo "genome:          			$genome"
	echo "reads:           			$read1"
	echo "                 			$read2"
	echo "mapped illumina:			$bamOut"
	echo "new genome (prefix):		$updatedGenome"
	echo "iteration:				$iteration"

	### Pilon requires mapping short reads to the genome 
	### Pilon then uses this mapping to infer sites of confidence, or lackthereof, to update the genome
	### As with every Broad software, Pilon requires an index for mappings

	bwa index $genome
	bwa mem -t 24 $genome $read1 $read2 > $samOut
	echo "Finished mapping Illumina reads to genome"

	samtools view -Suh $samOut | samtools sort -@ 24 -O BAM -o $bamOut -
	echo "Finished sorting bam file"

	samtools index -b $bamOut
	echo "Finished making index of bam file"

	echo "Starting pilon now"
	module load java/1.8.0_121
	java -Xmx60G -jar /global/scratch/dmai/projects/nasuta_group/albom_mini_assembly/pilon/pilon-1.22.jar --threads 24 --genome $genome --frags $bamFile --output test --outdir $genomeDir --output $updatedGenome

	echo "Pilon run finished"
	echo ""; echo ""
}
export -f runPilon

## runRacon requires 5 positional arguments, each separated by ONE(1) space
## 1: reference directory (do not include trailing /)			## 2: reference genome (STRING, do not include directory) 
## 3: read directory (do not include last /) 					## 4: long reads (STRING, do not include directory)  
## 5: racon iteration (INT)

### RACON for Dpse124Y 25X repeats normal ###
#runRacon /global/scratch/alisonn/04_asm/canu_124Y_repeats_feb2019 dpse124Y_canu_repeats_25x.contigs.fasta /global/scratch/alisonn/00_rawData/nanopore_v2 run15_run25_run27_all_pass_nanopore_het_Y_unmapped_singleLine.fastq.gz 1 
#runRacon /global/scratch/alisonn/04_asm/canu_124Y_repeats_feb2019 dpse124Y_canu_repeats_25x.racon_1.contigs.fasta /global/scratch/alisonn/00_rawData/nanopore_v2 run15_run25_run27_all_pass_nanopore_het_Y_unmapped_singleLine.fastq.gz 2 
#runRacon /global/scratch/alisonn/04_asm/canu_124Y_repeats_feb2019 dpse124Y_canu_repeats_25x.racon_2.contigs.fasta /global/scratch/alisonn/00_rawData/nanopore_v2 run15_run25_run27_all_pass_nanopore_het_Y_unmapped_singleLine.fastq.gz 3

## runPilon requires 4 positional arguments separated by space
## 1: genome directory 								2: genome (STRING, do not include directory)   				3: read directory 
## 4: short read pair 1 (STR, do not include dir)   5: short read pair 2 (same as 2) 					  		6: pilon iteration

#runPilon /global/scratch/alisonn/04_asm/canu_124Y_repeats_feb2019 dpse124Y_canu_repeats_25x.racon_3.contigs.fasta /global/scratch/alisonn/00_rawData/dpse-gDNA DBCC035B8_S64_L008_R1_001_val_1.fq.gz DBCC035B8_S64_L008_R2_001_val_2.fq.gz 1 

