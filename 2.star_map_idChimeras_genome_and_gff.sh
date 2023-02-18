#!/bin/bash

# v4 AP 26 Feb 2021

#dependencies
#STAR needs to be installed
#need to have base_genomics env active
#modify script if SE reads

#usage
#sbatch --account=palab -N 1 star_map.sh read1 read2 STAR_genome #for HPC
 
# e.g.
# sbatch --account=palab -N 1 star_map.sh S01_WT_Empty_Mock_1.fq.gz S01_WT_Empty_Mock_2.fq.gz STAR_genome #for HPC


##### read in command line
read1=$1
read2=$2
STAR_genome=$3
gtf=$4

prefix=$(ls $read1 | sed 's/_1_val_1.fq.gz//g' | sed 's/_2_val_2.fq.gz//g')

STAR --readFilesIn $read1 $read2 --outFileNamePrefix $prefix --genomeDir $STAR_genome --sjdbGTFfile $gtf --readFilesCommand zcat --runThreadN 12 --limitBAMsortRAM 1124320000 --chimScoreJunctionNonGTAG -1 --outSAMtype BAM SortedByCoordinate --chimOutType Junctions SeparateSAMold --alignSJDBoverhangMin 5 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 5 --outMultimapperOrder Random --outSAMattributes NH HI AS nM NM XS ch --chimSegmentMin 10 --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreSeparation 7 --chimSegmentReadGapMax 3 --chimFilter None --twopassMode None --alignSJstitchMismatchNmax 5 -1 5 5 --chimMainSegmentMultNmax 10  --outSAMattrRGline ID:my-read-group