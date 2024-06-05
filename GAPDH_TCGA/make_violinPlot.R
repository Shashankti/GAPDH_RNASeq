#GDC portal method

library(TCGAbiolinks)
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
library(ggplot2)
library(cowplot)
library(viridis)
library(tidyverse)
library(rstatix)
library(ggrepel)


TCGAbiolinks:::getProjectSummary("TCGA-SKCM")
options(timeout = '40000')



query_TCGA <- GDCquery(project = "TCGA-SKCM",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts"
)

GDCdownload(query_TCGA, method = "api", files.per.chunk = 100,
            directory = "~/TCGA_Data/")
rna_df <- GDCprepare(query_TCGA, directory = "~/TCGA_Data/")
rna_df <- assay(rna_df)

rownames(rna_df) <- sub('\\.[0-9]*$', '', rownames(rna_df))

other_merged_counts <- na.omit(merge(rna_df, merged_counts, by="row.names", all = TRUE))
other_merged_counts <- other_merged_counts %>% remove_rownames() %>% column_to_rownames(var = "Row.names")
other_merged_counts <- other_merged_counts[,which(colnames(other_merged_counts) %in% colnames(file2))]
other_merged_counts <- other_merged_counts[,-(192:217)]

col_data_2 = col_data
rownames(col_data_2)[1:150] <- colnames(other_merged_counts)[1:150] 


# Step 2: construct a DESeqDataSet object
dds_other <- DESeqDataSetFromMatrix (countData = round(other_merged_counts), colData = col_data_2 , design = ~ Type + Source)

# keep_2f <- rowSums(counts(dds_2f)) >= 10
# dds_2f <- dds_2f [keep_2f,]
# dds_2f

#get normalised counts information
dds_2f <- estimateSizeFactors(dds_2f)
sizeFactors(dds_2f)
dds_other$Type <- relevel(dds_other$Type,ref="Control")

# Step 3: Run DESeq ----------------------

dds_other <- DESeq(dds_other)
resultnames_2f <- resultsNames(dds_2f)
resultnames_2f

result_other_TCGA <- results(dds_other, contrast = c("Type", "Tumor", "Control")) 
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(result_TCGA), mart= mart)
result_other_TCGA$GeneName <- gene_IDs[match(rownames(result_other_TCGA), gene_IDs[,1]),2]
result_other_TCGA <- as.data.frame(result_other_TCGA)
result_other_TCGA$log2FoldChange <- as.numeric(result_other_TCGA$log2FoldChange)

write.table(result_TCGA, file = "All_DEG.tsv", quote = F, sep = "\t", row.names = T)
#read in non canonical genes
non_canonical_list <- read.table("NonCanonical.txt", header = F)

#plots
plot3 <- EnhancedVolcano::EnhancedVolcano(result_other_TCGA[result_other_TCGA$GeneName %in% non_canonical_list$V1,],
                                          lab = result_other_TCGA[result_other_TCGA$GeneName %in% non_canonical_list$V1,]$GeneName,
                                          x = 'log2FoldChange',
                                          y = 'padj',
                                          selectLab = c("GAPDH", "RPL18A","RPS10","RPL15","MRPL32"),
                                          title = 'Subsetting NonCanonical RBPs',
                                          subtitle = 'Differential Expressions',
                                          caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                                          pCutoff = 0.05,
                                          FCcutoff = 0.5,
                                          pointSize = 1.0,
                                          labSize = 3.0,
                                          colAlpha = 1,
                                          legendPosition = 'right',
                                          legendLabSize = 6,
                                          legendIconSize = 1.0,
                                          drawConnectors = TRUE,
                                          widthConnectors = 0.75,
                                          labCol = 'black',
                                          labFace = 'bold',
                                          boxedLabels = TRUE,
                                          colConnectors = 'black')
tiff("NoncanonicalRBPs_Subset.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(plot1)
dev.off()
png("NoncanonicalRBPs_Subset.png",width = 35,height = 21,units = 'cm',res = 300)
plot(plot1)
dev.off()


#limma pipeline
Type <- as.factor(col_data$Type)
Source <- as.factor(col_data$Source)
design = model.matrix(~Type+Source)
head(design)


dge = DGEList( # creating a DGEList object
  counts=final_merged_counts,
  samples=col_data,
  genes=result_TCGA$GeneName)

keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) 


dge = calcNormFactors(dge,method="TMM")
v = voom(dge,design,plot=TRUE)

fit = lmFit(v, design)
fit = eBayes(fit)


topGenes = topTable(fit, coef=1, sort.by="p", n=Inf)
print(topGenes)

limma_res <- list(voomObj = v, 
                  fit = fit,
                  topGenes = topGenes)

plot_PCA = function(voomObj, condition_variable){
  group = factor(voomObj$targets[, condition_variable])
  pca = prcomp(t(voomObj$E))
  # Take PC1 and PC2 for the plot
  plot(pca$x[,1:2],col=group, pch=19)
  # include a legend for points
  legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
  return(pca)
}

res_pca = plot_PCA(limma_res$voomObj, "Type")

#make violin plot
topGenes$adj.P.Val <- as.numeric(topGenes$adj.P.Val)

#plots
limma_vol <- EnhancedVolcano::EnhancedVolcano(topGenes[topGenes$genes %in% non_canonical_list$V1,],
                                          lab = topGenes[topGenes$genes %in% non_canonical_list$V1,]$genes,
                                          x = 'logFC',
                                          y = 'adj.P.Val',
                                          selectLab = c("GAPDH", "RPL18A","RPS10","RPL15","MRPL32"),
                                          title = 'Subsetting NonCanonical RBPs',
                                          subtitle = 'Differential Expressions',
                                          caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                                          pCutoff = 0.05,
                                          FCcutoff = 0.5,
                                          pointSize = 1.0,
                                          labSize = 3.0,
                                          colAlpha = 1,
                                          legendPosition = 'right',
                                          legendLabSize = 6,
                                          legendIconSize = 1.0,
                                          drawConnectors = TRUE,
                                          widthConnectors = 0.75,
                                          labCol = 'black',
                                          labFace = 'bold',
                                          boxedLabels = TRUE,
                                          colConnectors = 'black')

EnhancedVolcano::EnhancedVolcano(topGenes,
                                 lab = topGenes$genes,
                                 x = 'logFC',
                                 y = 'adj.P.Val',
                                 selectLab = c("GAPDH", "RPL18A","RPS10","RPL15","MRPL32"),
                                 title = 'Subsetting NonCanonical RBPs',
                                 subtitle = 'Differential Expressions',
                                 caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                                 pCutoff = 0.05,
                                 FCcutoff = 0.5,
                                 pointSize = 1.0,
                                 labSize = 3.0,
                                 colAlpha = 1,
                                 legendPosition = 'right',
                                 legendLabSize = 6,
                                 legendIconSize = 1.0,
                                 drawConnectors = TRUE,
                                 widthConnectors = 0.75,
                                 labCol = 'black',
                                 labFace = 'bold',
                                 boxedLabels = TRUE,
                                 colConnectors = 'black')




# Transpose and make it into a matrix object
d_mat = as.matrix(t(limma_res$voomObj$E))





clinical = data@colData

dim(clinical)

expr_diseased = d_mat[Tumor_samples,gene_id]
expr_healthy = d_mat[Control_samples, gene_id]
boxplot(expr_diseased, expr_healthy,
        names=c("Diseased", "Healthy"), main="Distribution of gene expression")

clin_df = clinical[,c("patient",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "gender",
                     "tumor_stage")]

#make violin plot
gene_id = "ENSG00000111640"
gene_name = "GAPDH"
Control_samples = colnames(GAPDH_all[,which(colnames(normalized_counts) %in% rownames(col_data[col_data$Type=="Control",]))])
Tumor_samples = colnames(GAPDH_all[,-which(colnames(normalized_counts) %in% rownames(col_data[col_data$Type=="Control",]))])
GAPDH_all[,which(colnames(normalized_counts) %in% rownames(col_data[col_data$Type=="Control",]))]
raw_GAPH <- t(final_merged_counts[rownames(final_merged_counts) == gene_id,])


GAPDH_all <- as.data.frame(GAPDH_all)
GAPDH_all$Type <- ifelse(rownames(GAPDH_all) %in% Control_samples, "Normal", "AML")
GAPDH_all$Raw <- raw_GAPH
colnames(GAPDH_all) <- c("Normalized", "Type","Raw")

p1 <- ggplot(GAPDH_all, aes(x = fct_reorder(Type, Raw), y = Raw, fill = Type))+
  geom_violin(trim = FALSE, width=0.4)+
  theme_cowplot()+
  geom_boxplot(width = 0.1, color = "grey", alpha=0.5)+
  theme(legend.position = "none",
        plot.title = element_text(size = 11))+
  ggtitle("GAPDH Raw counts")+
  geom_jitter(height = 0, width = 0.1,alpha=0.7)+
  scale_fill_viridis_d()+
  ylab("Raw Counts") +
  xlab("Type")

p1 <- p1+
  annotate("text",x = 0.58, xend = 1.9, y= 5.2, label = "LFC = 1.353397")+
  annotate("text",x = 0.58, xend = 1.9, y= 5.1, label = "pvalue =  4.98e-08")
p1 <- p1+
  annotate("text",x = 0.58, xend = 1.9, y= 200000, label = "LFC = 1.353397")+
  annotate("text",x = 0.58, xend = 1.9, y= 190000, label = "pvalue = 4.98e-08")
png("Gapdh_raw.png",width = 35,height = 21,units = 'cm',res = 300)
plot(p1)
dev.off()

#wilcoxons rank sum test
GAPDH_all$Log_Norm <- NULL
GAPDH_all$Log_Raw <- NULL
stat.test <- GAPDH_all %>%
  wilcox_test(Normalized ~ Type) %>%
  add_significance()
stat.test

GAPDH_all %>% wilcox_effsize(Raw ~ Type)
stat.test <- stat.test%>%add_xy_position(x="Type")

fold.change <- log(mean(GAPDH_all[GAPDH_all$Type=="AML",3])/mean(GAPDH_all[GAPDH_all$Type=="Normal",3]),base=2)



  annotate("text",x=-1.8,y=0.97,label = "(FC=0.005864923)", color = "red") +
  annotate("text",x=-1.8,y=0.90, label = "324 Gal Genes",color = "#7dad02") + 
  annotate("text",x=-1.8,y=0.87, label = "(FC=-0.2044665,", color = "#7dad02") + 
  annotate("text",x=-1.8,y=0.84, label = "P = 6.44e-04)", color = "#7dad02") + 
  annotate("text",x=-1.8,y=0.78, label = "449 Low Genes", color = "#c57ef8") +
  annotate("text",x=-1.8,y=0.75,label = "(FC=-0.1975091,", color = "#c57ef8") +
  annotate("text",x=-1.8,y=0.72,label = "P=9.35e-06)",color = "#c57ef8") +
  annotate("text",x=-1.8,y=0.67,label= "375 Norm Genes", color = "#00bfc4") +
  annotate("text",x=-1.8,y=0.64,label = "(FC=-0.2585118,", color = "#00bfc4") + 
  annotate("text",x=-1.8,y=0.61,label = "P=5.043-06", color = "#00bfc4")


