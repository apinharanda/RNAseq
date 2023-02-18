#!/bin/sh 

#####feb14
#v1

##this is a script to count features in bam for RNAseq analysis
#Read the documentation of htseqcount
#pay special attention to stranded option!! if RNAseq data is stranded usually the option is "reverse" not a much more intuitive "stranded"
#can change the feature to exon or whatever
##IMPORTANT: below the option is no for stranded, bc yes, the default is v rarely correct (it will be either no, or reverse most of the times)

#have base_genomics activated

#can be run like 
#sbatch --account=palab --time=05:00:00 htseqcount.sh example.bam example_counts.txt

#make bam index before counting - samtools index bam

###read from command line

name=$1 #bam name
gtf=$2 #gtf with complete path, if you only have a gff you need to convert it to a gtf, use gffread my.gff3 -T -o my.gtf
outname=$3


htseq-count --format=bam --type=gene --order=pos --stranded=no --idattr=gene $name  $gff > $outname
