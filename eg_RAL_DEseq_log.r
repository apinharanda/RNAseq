library("DESeq2")
library("ggplot2")
library("ggforce")
library("RColorBrewer")
library("pheatmap")

#R --vanilla --slave '--args sample_table_info_.txt normalised_counts_.txt rld_ddsHTseq_sampleDists_Matrix_.pdf genes_PCs_.txt PCA_ercentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix_.pdf heat_mat_rld_ddsHTseq_geneCluster_geneName_.pdf mat_rld_ddsHTseq_TOPgenes_.txt mat_rld_ddsHTseq_categories_TOPgenes_.txt results_contrast_.txt' <DESeq2_move.R > DESeq2_move_.log

#Rscript DESeq2_move.R sample_table_info_.txt normalised_counts_.txt rld_ddsHTseq_sampleDists_Matrix_.pdf genes_PCs_.txt PCA_ercentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix_.pdf heat_mat_rld_ddsHTseq_geneCluster_geneName_.pdf mat_rld_ddsHTseq_TOPgenes_.txt mat_rld_ddsHTseq_categories_TOPgenes_.txt results_contrast_.txt  > DESeq2_move_.log

rm(list = ls())

args <- commandArgs()

cat(args[5],args[6],args[7],args[8],args[9],args[10],args[11],args[12],args[13],"\n")

sample_table_info_out = (args[5]) #sample_table_info.txt
normalised_counts_out = (args[6]) #normalised_counts.txt
rld_ddsHTseq_sampleDists_Matrix_out = (args[7]) #rld_ddsHTseq_sampleDists_Matrix.pdf
genes_PCs_out = (args[8]) #genes_PCs.txt
PCA_ercentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix_out = (args[9]) #PCA_ercentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix.pdf
heat_mat_rld_ddsHTseq_geneCluster_geneName_out = (args[10]) #heat_mat_rld_ddsHTseq_geneCluster_geneName.pdf
mat_rld_ddsHTseq_TOPgenes_out = (args[11]) #mat_rld_ddsHTseq_TOPgenes.txt
mat_rld_ddsHTseq_categories_TOPgenes_out = (args[12]) #mat_rld_ddsHTseq_categories_TOPgenes.txt
results_contrast_out = (args[13]) #results_contrast.txt



directory <- "/moto/palab/users/apinharanda/caudalHorn/data_firefly/raw_reads"

sampleTable <- read.table (sample_table_info, header = TRUE)
#sampleTable <- read.table ("sampleTable.txt", header = TRUE)


ddsHTseq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = directory,  
                                      design = ~Species:Replicate + Species*Tissue)

head (sampleTable)


#get normalised counts


dds <- estimateSizeFactors(ddsHTseq)
normalised_counts <- counts(dds, normalized=TRUE)

write.table(as.data.frame(normalised_counts), file=normalised_counts_out,sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
#write.table(as.data.frame(normalised_counts), file="normalised_counts_out.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)

##estimate dispersion
dds_disp <- estimateDispersions(dds)
head(dispersions(dds_disp))



rld_ddsHTseq <- rlog(ddsHTseq,  blind=FALSE)


#genes with low counts seem to be excessively variable on the ordinary logarithmic scale
#while the rlog transform compresses differences for genes for which the data cannot provide good information anyway

#sample distances
#access similarity between samples
#dist to calculate Eucledian distance between samples
#use rlog data to avoid that the distance measure is dominated by a few highly variable genes and has a roughly equal contribution form all genes


rld_ddsHTseq_sampleDists <- dist( t( assay(rld_ddsHTseq) ) )

rld_ddsHTseq_sampleDists_Matrix <- as.matrix( rld_ddsHTseq_sampleDists )

rownames(rld_ddsHTseq_sampleDists_Matrix) <- paste( rld_ddsHTseq$Species, rld_ddsHTseq$Tissue, rld_ddsHTseq$Replicate, sep="-" )

colnames(rld_ddsHTseq_sampleDists_Matrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf (rld_ddsHTseq_sampleDists_Matrix_out)
#pdf ("rld_ddsHTseq_sampleDists_Matrix_out.pdf")
pheatmap(rld_ddsHTseq_sampleDists_Matrix,
         clustering_distance_rows=rld_ddsHTseq_sampleDists,
         clustering_distance_cols=rld_ddsHTseq_sampleDists,
         col=colors)
dev.off()


###proportion of the variance

rv = rowVars(assay(rld_ddsHTseq)) 
select = order(rv, decreasing=TRUE)[seq_len(min(500, length(rv)))]
pc = prcomp(t(assay(rld_ddsHTseq)[select,]))
loadings = as.data.frame(pc$rotation)
aload = abs(loadings)
genes_PCs <- sweep(aload, 2, colSums(aload), "/")

write.table(as.data.frame(genes_PCs), file=genes_PCs_out,sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
#
##PCA

pcaData_rld_ddsHTseq_sampleDists_Matrix <- plotPCA(rld_ddsHTseq, intgroup = c("Species","Tissue"), returnData = TRUE)
pcaData_rld_ddsHTseq_sampleDists_Matrix

percentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix <- round(100 * attr(pcaData_rld_ddsHTseq_sampleDists_Matrix, "percentVar"))



pdf(PCA_ercentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix_out)
#pdf("PCA_ercentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix_out.pdf")

ggplot(pcaData_rld_ddsHTseq_sampleDists_Matrix, aes(x = PC1, y = PC2, color = Species, shape=Tissue)) +
geom_point(size=2) +
theme_bw()+
xlab(paste0("PC1: ", percentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix[1], "% variance")) +
ylab(paste0("PC2: ", percentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix[2], "% variance")) +
theme(axis.line = element_line(colour = "black"), panel.border=element_rect(colour="black", size=1),aspect.ratio=0.6)

dev.off()



pdf(PCA_ercentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix_out)
#pdf("PCA_ercentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix_out.pdf")

ggplot(pcaData_rld_ddsHTseq_sampleDists_Matrix, aes(x = PC1, y = PC2, shape = Tissue,colour=Species)) +
geom_point(size=2) +
theme_bw()+
xlab(paste0("PC1: ", percentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix[1], "% variance")) +
ylab(paste0("PC2: ", percentVar_pcaData_rld_ddsHTseq_sampleDists_Matrix[2], "% variance")) +
theme(axis.line = element_line(colour = "black"), panel.border=element_rect(colour="black", size=1),aspect.ratio=0.6)

dev.off()




#######Gene clustering
#these are the top 20 genes
#that is enough to show in a paper and a decent number for you  to be able to manually see what they may be in flybase
#you can change the number to be however many you want

topVarGenes <- head(order(-rowVars(assay(rld_ddsHTseq))),50) #here change 50, for the number you want


mat_rld_ddsHTseq <- assay(rld_ddsHTseq)[ topVarGenes, ]
mat_rld_ddsHTseq <- mat_rld_ddsHTseq - rowMeans(mat_rld_ddsHTseq)
df_mat_rld_ddsHTseq <- as.data.frame(colData(rld_ddsHTseq)[,c("Species","Tissue")])

##this was done for the caudal horn paper where i wanted to include specific genes that overlaped QTL
#does not need to be done for a regular RNAseq analysis
#genes <- c("KWMTBOMO01772","KWMTBOMO02806","KWMTBOMO02816","KWMTBOMO02819","KWMTBOMO02823","KWMTBOMO02831","KWMTBOMO02843","KWMTBOMO02872","KWMTBOMO02879","KWMTBOMO02896","KWMTBOMO02897","KWMTBOMO02906","KWMTBOMO02915","KWMTBOMO02921","KWMTBOMO02923","KWMTBOMO02933","KWMTBOMO02944","KWMTBOMO02953","KWMTBOMO02955","KWMTBOMO02956","KWMTBOMO02957","KWMTBOMO02958","KWMTBOMO02960","KWMTBOMO02962","KWMTBOMO02979","KWMTBOMO02996","KWMTBOMO02997","KWMTBOMO03010","KWMTBOMO03030","KWMTBOMO03045","KWMTBOMO03047","KWMTBOMO04838","KWMTBOMO04865","KWMTBOMO04867","KWMTBOMO04884","KWMTBOMO04885","KWMTBOMO04895","KWMTBOMO04898","KWMTBOMO04899","KWMTBOMO04903","KWMTBOMO04904","KWMTBOMO04906","KWMTBOMO04907","KWMTBOMO04913","KWMTBOMO04914","KWMTBOMO04956","KWMTBOMO04965","KWMTBOMO04967","KWMTBOMO04973","KWMTBOMO04974","KWMTBOMO04985","KWMTBOMO04986","KWMTBOMO04990","KWMTBOMO05026","KWMTBOMO05498","KWMTBOMO05521","KWMTBOMO05535","KWMTBOMO05541","KWMTBOMO05547","KWMTBOMO05551","KWMTBOMO05559","KWMTBOMO05567","KWMTBOMO05587","KWMTBOMO05616","KWMTBOMO05655","KWMTBOMO05656","KWMTBOMO05692","KWMTBOMO05726","KWMTBOMO05734","KWMTBOMO05735","KWMTBOMO05737","KWMTBOMO05757","KWMTBOMO05758","KWMTBOMO05762","KWMTBOMO05764","KWMTBOMO05771","KWMTBOMO05805","KWMTBOMO05817","KWMTBOMO05830","KWMTBOMO05834","KWMTBOMO05861","KWMTBOMO05865","KWMTBOMO05871","KWMTBOMO05883","KWMTBOMO05913","KWMTBOMO05940","KWMTBOMO05941","KWMTBOMO05944","KWMTBOMO05947","KWMTBOMO05949","KWMTBOMO05957","KWMTBOMO05959","KWMTBOMO05983","KWMTBOMO05994","KWMTBOMO06011","KWMTBOMO06018","KWMTBOMO06031","KWMTBOMO06064","KWMTBOMO12587","KWMTBOMO12625","KWMTBOMO12649","KWMTBOMO12650","KWMTBOMO12651","KWMTBOMO12662","KWMTBOMO12686","KWMTBOMO12689","KWMTBOMO12694","KWMTBOMO12695","KWMTBOMO12696","KWMTBOMO12705","KWMTBOMO12709","KWMTBOMO12711","KWMTBOMO15632","KWMTBOMO15652","KWMTBOMO15653","KWMTBOMO15660","KWMTBOMO15707","KWMTBOMO15720","KWMTBOMO15721","KWMTBOMO15722","KWMTBOMO15726")

#mat_rld_ddsHTseq <- assay(rld_ddsHTseq)[genes,]
#mat_rld_ddsHTseq <- mat_rld_ddsHTseq - rowMeans(mat_rld_ddsHTseq)
#df_mat_rld_ddsHTseq <- as.data.frame(colData(rld_ddsHTseq)[,c("Tissue","Species")])
#

#mat_121_rld_ddsHTseq_mean <- read.table ("mat_121_rld_ddsHTseq.txt")
#or
#mat_121_rld_ddsHTseq_mean<- read.table ("normalised_counts_121.txt")
#or
#mat_121_rld_ddsHTseq_mean<- read.table ("normalised_counts_121_log2.txt")
#df <- mat_121_rld_ddsHTseq_mean[1,]

#df_selected <- read.table ("df_selected.txt")

#pdf ("heatmap_overlpsQTL_121_final.pdf")
#pheatmap(mat_121_rld_ddsHTseq_mean, annotation_col=df_selected,fontsize = 4, cluster_rows=TRUE, cluster_cols=FALSE,fontsize_row = 4, fontsize_col = 5,color=colorRampPalette(c("blue","white","red"))(200))
#dev.off()

#pdf ("heatmap_overlpsQTL_121_wnt6_final.pdf")
#pheatmap(df, annotation_col=df_selected,fontsize = 5, fontsize_row = 4, cluster_rows=FALSE, cluster_cols=FALSE,fontsize_col = 5,color=colorRampPalette(c("blue","white","red"))(200))
#dev.off()

#KWMTBOMO01772 7.835006       6.687298    6.811473       6.804376
#228.3346      103.05695   112.32014       111.7690
#write.table(as.data.frame(mat_rld_ddsHTseq), file="mat_rld_ddsHTseq.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)



pdf (heat_mat_rld_ddsHTseq_geneCluster_geneName_out)
#pdf ("heat_mat_rld_ddsHTseq_geneCluster_geneName_out.pdf")
pheatmap(mat_rld_ddsHTseq, annotation_col=df_mat_rld_ddsHTseq,fontsize = 6, fontsize_row = 5, fontsize_col = 7)
dev.off()


mat_rld_ddsHTseq_df <- as.data.frame(mat_rld_ddsHTseq)
write.table(as.data.frame(mat_rld_ddsHTseq_df), file=mat_rld_ddsHTseq_TOPgenes_out,sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
write.table(as.data.frame(df_mat_rld_ddsHTseq), file=mat_rld_ddsHTseq_categories_TOPgenes_out,sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)

#write.table(as.data.frame(mat_rld_ddsHTseq_df), file="mat_rld_ddsHTseq_TOPgenes_out.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
#write.table(as.data.frame(df_mat_rld_ddsHTseq), file="mat_rld_ddsHTseq_categories_TOPgenes_out.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)



##########PART 2 
#DIFFERENTIAL EXPRESSION ANALYSIS

dds_DESeq_ddsHTseq <- DESeq(ddsHTseq)

#this makes low the baseline
#so positivee fold change are high group being more expressed than low group

ddsHTseq$Tissue <- relevel(ddsHTseq$Tissue, "NonCaudal")
ddsHTseq$Tissue <- relevel(ddsHTseq$Tissue, "Caudal")

results_contrast <- results (dds_DESeq_ddsHTseq, contrast = c("Tissue", "Caudal", "NonCaudal"))

#this is the table with all the results
#the correct pvalue to look at is the adjusted pvalue (adpval), not the pval
#the adjusted is corrected for multiple testing
write.table(as.data.frame(results_contrast), file="results_contrast_out_caudal_vs_nonCaudal.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)



#this makes low the baseline
#so positivee fold change are high group being more expressed than low group

ddsHTseq$Species <- relevel(ddsHTseq$Species, "Bmori")

results_contrast <- results (dds_DESeq_ddsHTseq, contrast = c("Species", "Bmori", "Bmand"))

#this is the table with all the results
#the correct pvalue to look at is the adjusted pvalue (adpval), not the pval
#the adjusted is corrected for multiple testing
write.table(as.data.frame(results_contrast), file="results_contrast_out_Bmori_vs_Bmand.txt",sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)

