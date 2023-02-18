#!/bin/bash
#Feb 21 2022
#can be run in the head node
#activate base_genomics
#this is just an example


STAR --runMode genomeGenerate --genomeDir BomoRef --genomeFastaFiles ../../reference_genomes/Bomo_genome_assembly.fa --sjdbGTFfile ../../reference_genomes/Bomo_gene_models.gff3 --sjdbOverhang 99 --outFileNamePrefix BomoSTAR --genomeSAindexNbases 13