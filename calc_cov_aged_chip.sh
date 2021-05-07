#!/bin/bash
# Job name:
#SBATCH --job-name=bedtl
#SBATCH --account=co_bachtrog
#SBATCH --partition=savio2
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
# general use for getting read coverage

# Stop the whole script on CTRL+C
int_handler()
{
    echo "Interrupted."
    # Kill the parent process of the script.
    kill $PPID
    exit 1
}
trap 'int_handler' INT

export PATH=$PATH:/global/home/groups/co_bachtrog/programs/samtools-1.5
export PATH=$PATH:/global/home/groups/co_bachtrog/programs/bedtools2-2.19.1/bin
runCounts () {

    sample=$(echo $1 | cut -f 1 -d '#')
    sample_name=$(echo $sample | cut -f 1 -d '.')
    sample_type=$(echo $sample | cut -f 2 -d '.')
    size=100 ## any resolution you want

    bamDir="/global/scratch/alisonn/02_mapped/dmir-aging/chip"
    covDir="/global/scratch/alisonn/03_cov/dmir-aging/revisions/100bp"
    windows="/global/home/users/alisonn/01_ref/"
    sample_windows=/global/home/users/alisonn/01_ref/dmir/FinalApril2018_pilon.${size}bp_windows.bed
    spike_windows=/global/home/users/alisonn/01_ref/dmel/dmel-all-chromosome-r6.14.${size}bp_windows.bed

    inBam=${bamDir}/${sample}
    outBed=${covDir}/${sample_name}.${sample_type}.${size}bp.temp.counts
    sortedBed=${covDir}/${sample_name}.${sample_type}.${size}bp.counts

    if [ "$sample_type" = "mappedtodmel" ] ; then
        windows=$spike_windows
    else
        windows=$sample_windows
    fi

    echo "Making counts for this bam file: $inBam "
    echo "Using this window: $windows"
    echo "The sorted file will be named: $sortedBed"
	echo "This is the sample type: $sample_type"
    echo ""
	echo "" 

    #  compute the counts for each bam file
    bedtools coverage -counts -abam $inBam -b $windows > $outBed
    sort -k1,1 -k2,2n $outBed > $sortedBed
    rm $outBed

}

export -f runCounts
## example usage of runCounts()
#/global/home/users/alisonn/parallel runCounts ::: MSH22_F_5d_H3K9me3_rep1_ChIP.mappedtodmel.sorted.bam \
#MSH22_F_5d_H3K9me3_rep1_ChIP.mappedtodmir.sorted.bam \
#MSH22_F_5d_H3K9me3_rep1_input.mappedtodmel.sorted.bam \
#MSH22_F_5d_H3K9me3_rep1_input.mappedtodmir.sorted.bam \
#MSH22_F_5d_H3K9me3_rep2_ChIP.mappedtodmel.sorted.bam \
#MSH22_F_5d_H3K9me3_rep2_ChIP.mappedtodmir.sorted.bam \
#MSH22_F_5d_H3K9me3_rep2_input.mappedtodmel.sorted.bam \
#MSH22_F_5d_H3K9me3_rep2_input.mappedtodmir.sorted.bam \
#MSH22_F_80d_H3K9me3_rep3_ChIP.mappedtodmel.sorted.bam 

