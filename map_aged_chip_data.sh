#!/bin/bash
# Job name:
#SBATCH --job-name=bt2
#SBATCH --account=co_bachtrog
#SBATCH --partition=savio2_bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24 
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

# mapping reads to both sample and spike reference genomes with bowtie2 parameters from Brown et al 2020

# Stop the whole script on CTRL+C
int_handler() {
    echo "Interrupted."
    # Kill the parent process of the script.
    kill $PPID
    exit 1
}
trap 'int_handler' INT

runMapping () {
    readsDir="/global/scratch/alisonn/00_rawData/dmir-aging/subsampled"
    outDir="/global/scratch/alisonn/02_mapped/dmir-aging/revisions"
    refDir="/global/home/users/alisonn/01_ref"

    ### SET UP THE SAMPLES AND GENOMES FOR MAPPING
    input=$1
    sampleL=$( echo $1 | cut -f 1 -d \# )
    sampleR=$( echo $1 | cut -f 2 -d \# )
    prefix=$( echo $1 | cut -f 3 -d \# )
	sex=$( echo $prefix | cut -f 2 -d '_' )

    ### SET UP MAPPING PROPERLY
    sampleType=$( echo $1 | cut -f 4 -d \# )
    referenceGenome=""
    spikeReference=dmel-all-chromosome-r6.14
    maleReference=FinalApril2018_pilon
	femaleReference=FinalApril2018_pilon_femaleGenome
	
    if [ "$sampleType" = "mappedtodmel" ] ; then
        referenceGenome=$spikeReference
    elif [ "$sex" == "F" ] ; then
        referenceGenome=$femaleReference
	else
		referenceGenome=$maleReference
    fi

    fastq1=${readsDir}/${sampleL}
    fastq2=${readsDir}/${sampleR}
    mapSam=${outDir}/${prefix}.${sampleType}.sam
    outBam=${outDir}/${prefix}.${sampleType}.sorted.bam

    echo "The command line input was this $input"
    echo "Which is parsed as this: $fastq1"
    echo "$fastq2"
    echo "$outBam"
    echo "Mapping to this $referenceGenome for the $sampleType"
    echo " "
    echo " "
    echo " "

    bowtie2 -p 8 -D 15 -R 2 -N 0 -L 22 -i S,1,0.50 --no-unal --no-1mm-upfront -x $referenceGenome -1 ${fastq1} -2 ${fastq2} -S ${mapSam}
    samtools view -Suh $mapSam | samtools sort -T ${prefix} -O bam -@ 8 -o ${outBam} -

    rm $mapSam
}

export -f runMapping

## Make index files for genomes
#bowtie2-build /global/home/users/alisonn/01_ref/dmir/FinalApril2018_pilon.fasta FinalApril2018_pilon
#bowtie2-build /global/home/users/alisonn/01_ref/dmir/FinalApril2018_pilon_femaleGenome.fasta FinalApril2018_pilon_femaleGenome
#bowtie2-build /global/home/users/alisonn/01_ref/dmir/dmir-all-chromosome-r6.14.fasta dmir-all-chromosome-r6.14

## Example list of input format for runMapping() 
## can also input a txt file for GNU parallel
#/global/home/users/alisonn/parallel -j 3 runMapping ::: DBAN002A_S46_L008_R1_001_val_1.fq.gz#DBAN002A_S46_L008_R2_001_val_2.fq.gz#MSH22_F_5d_H3K9me3_rep1_ChIP#mappedtodmir \
#DBAN002B_S47_L008_R1_001_val_1.fq.gz#DBAN002B_S47_L008_R2_001_val_2.fq.gz#MSH22_M_5d_H3K9me3_rep1_ChIP#mappedtodmir \
#DBAN002C_S48_L008_R1_001_val_1.fq.gz#DBAN002C_S48_L008_R2_001_val_2.fq.gz#MSH22_F_98d_H3K9me3_rep1_ChIP#mappedtodmir \
#DBAN002D_S49_L008_R1_001_val_1.fq.gz#DBAN002D_S49_L008_R2_001_val_2.fq.gz#MSH22_M_98d_H3K9me3_rep1_ChIP#mappedtodmir \
#DBAN002G_S52_L008_R1_001_val_1.fq.gz#DBAN002G_S52_L008_R2_001_val_2.fq.gz#MSH22_F_5d_H3K9me3_rep1_input#mappedtodmir \
#DBAN002H_S53_L008_R1_001_val_1.fq.gz#DBAN002H_S53_L008_R2_001_val_2.fq.gz#MSH22_M_5d_H3K9me3_rep1_input#mappedtodmir 
