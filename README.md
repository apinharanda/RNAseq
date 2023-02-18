# RNAseq
Scripts and logs to do RNAseq analysis

Summary of RNAseq analysis pipeline
#############Setting up references

#create index on fasta to be used with samtools 

samtools faidx [fasta file]

#create a bwa index 

bwa index [fasta file]

#create the picard dictionary

java -jar picard.jar CreateSequenceDictionary R= [fasta file] O= [fasta file].dict 

#download samples (SRA or sequencing centre)


#############Working with the sequences

1. QC of the run
fastqc .gz 
then download .html and look at it

2. Read trimming
Trim adapters and reads that have too small or have too many errors
Using cutadapt & trimgalore

3. Fastqc again
Step 1 does not need to be necessarily done if this one is done

4. Make STAR genome dir

5. Map reads
With STAR reads to genome

e.g.
STAR --genomeDir ../genomeDir/ISO1_MT/ --sjdbGTFfile ../genomeDir/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --readFilesCommand zcat --outFileNamePrefix RAL21_fem_high_1 --readFilesIn RAL21_fem_high_1_trimmed.fq.gz

6. Sort bam using
samtools view

7. Check read depth and general stats
with samtools depth (but also samtools flagstat) or Picard CollectAlignmentSummaryMetrics

8. Count features (i.e. genes)
Using htseq

e.g.
htseq-count --format=bam --type=gene --order=pos --idattr=gene RAL21_fem_high_1Aligned.sortedByCoord.out.bam  ../genomeDir/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff > RAL21_fem_high_1_htseq_counts.txt

9. Differential expression analysis
With DESeq2


##Good papers & resources to learn how to design a matrix and generally understand how DESeq 2 works
http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-deseqdataset-object-sample-information-and-the-design-formula
http://rafalab.github.io/pages/harvardx.html
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/



