#!/bin/bash

# v5 AP 13 Feb 2022

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
prefix=$2
STAR_genome=$3

STAR --readFilesIn $read1  --outFileNamePrefix $prefix --genomeDir $STAR_genome  --readFilesCommand zcat --runThreadN 12 --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS nM NM XS ch --outSAMattrRGline ID:my-read-group