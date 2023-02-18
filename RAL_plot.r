###v5
#changed error on plotting transformed log2
#removed -/+ 1 for FC plotting 
#add other models
#changed y limits for plotting
#removed other models that were not necessary/wrong
#changed y axis out of faceting
#removed box plots

#2May 2021
#plot RAL normalised counts + fold change males & females

rm(list = ls())

#need to have the following installed
#tested on R 4.0.2

library("ggplot2")
library("tidyr") #v0.8.3
library("dplyr") #v0.8.3
library("ggpubr")

# to run this script
# R --vanilla --slave '--args Gene_name normalised_counts male_log2FC female_log2FC all_log2FC ind_log2FC' <RAL_plot.R; rm -rf Rplots.pdf

#need to input:
#1) gene name,
#2) normalised count file, normalised_counts_v2.txt
#3) male contrast, results_contrast_males.txt
#4) female contrast, results_contrast_females.txt

#for example
#R --vanilla --slave '--args corolla normalised_counts_v2.txt results_contrast_v2_ind.txt results_contrast_v2_cnd_plus_cnd_inter_group.txt' <RAL_plot.r; rm -rf Rplots.pdf

#assumes all files are in your path

args <- commandArgs()
cat(args[5],args[6],args[7],args[8],"\n")

Gene_name=args[5]

normalised_counts=args[6]

ind_log2FC=args[7]

group_inter_cnd_log2FC=args[8]


#colours
cbpal <- c("#0072B2", "#D55E00", "#56B4E9", "#009E73", "#000000", "#E69F00", "#F0E442", "#CC79A7")
names(cbpal)<-c("blue","red","lightblue","green","black","orange","yellow","magenta")

##plot gene counts
#change the normalised_counts_v2 to add Gene to the first column
#then it can be imported

normalised_counts <- read.table(normalised_counts,header=TRUE)

#do the same to the contrasts
#change the results_contrast_v2 to add Gene to the first column
#then it can also be imported

results_contrast_v2_ind<- read.table(ind_log2FC,header=TRUE)
results_contrast_group_inter_cnd <- read.table(group_inter_cnd_log2FC,header=TRUE)

results_contrast_v2_ind$sex<-"phenotype_sex_ind"
results_contrast_group_inter_cnd$sex <- "phenotype_sex"


results_contrast <- rbind (results_contrast_v2_ind,results_contrast_group_inter_cnd)


#Get normalised counts for gene of interest

gene <- normalised_counts[ which(normalised_counts$Gene==Gene_name), ]
gene_long <- pivot_longer(gene, cols = starts_with("RAL"), names_to = "Description", values_to = "Normalised_counts")
gene_long_sep <- separate(data = gene_long, col = Description, into = c("RAL", "sex","group","rep"), sep = "_")
gene_long_sep_df <- as.data.frame(gene_long_sep)


#Get adjusted pval

gene_pval <- results_contrast[ which(results_contrast$Gene==Gene_name), ]
gene_pval$FoldChange_plot <- 2 ^(abs(gene_pval$log2FoldChange)) 

ymax=max(gene_pval$FoldChange_plot)+0.5

#gene_pval$FoldChange_plot <- ifelse(gene_pval$log2FoldChange < 0, gene_pval$FoldChange*-1,ifelse(gene_pval$FoldChange >=0, gene_pval$FoldChange))
gene_pval$sign <- ifelse(gene_pval$log2FoldChange < 0, "Sensitive-HigherExpression",ifelse(gene_pval$FoldChange >=0, "Resistant-HigherExpression"))


##NOW THE ACTUAL PLOTS
#chage the gene and the titles etc

title_1=paste(Gene_name,"Normalized counts",sep="\n\n\n")

gene_long_sep_df
#title_1

gene_long_sep_df_plot=
ggplot(gene_long_sep_df) +
geom_point(aes(x=RAL, y=Normalised_counts, colour=sex),size=1.5,alpha=0.5) + 
facet_wrap(~group, scales="free_x")+
theme_bw() +
theme(aspect.ratio=1.3, legend.position="bottom")+
ggtitle(title_1)+
xlab("Strain")+
ylab("Nomalised counts")


gene_pval

gene_pval_plot=
ggplot(gene_pval,aes(y=FoldChange_plot, x=sex,colour=as.factor(padj)))+
geom_point(shape=95, size=3,stroke =15)+
facet_wrap(~sign, scales="free_x")+
theme_bw()+
theme(aspect.ratio=1.3,axis.text.x = element_text(angle=45,vjust=1,hjust=1))+
xlab("Model")+
ylab("FoldChange")+
ggtitle("Fold Change")+
labs(fill = "padj")+
geom_hline(aes(yintercept=1), colour="darkblue", linetype="dashed")+
scale_y_continuous(breaks=seq(1,ymax+0.5,1))+
scale_colour_manual(name="padj",values=c("#E69F00", "#009E73"))

chr_plot=ggarrange(gene_long_sep_df_plot,gene_pval_plot,nrow = 2)

title_2=paste(Gene_name,"Normalized_counts_logFC","pdf",sep=".")
#title_2

ggsave(filename=title_2, plot=chr_plot, dpi=500, units="cm", width=16.0, height=24.0)
