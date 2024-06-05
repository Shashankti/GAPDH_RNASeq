library(dplyr)
library(tidyverse)
library(cowplot)
library(ggplot2)
library(DESeq2)
library(biomaRt)
library(ggsignif)
library(PCAtools)
#read in old and new hsc counts file


old_hsc <- data.table::fread("old_hsc.tsv")
new_hsc <- data.table::fread("HSC_merged_counts_without_cas9.tsv")
combined_hsc <- na.omit(merge(new_hsc,old_hsc, by="ID",all=TRUE))
combined_hsc <- combined_hsc %>% remove_rownames() %>% column_to_rownames(var = "ID")
rownames(combined_hsc) <- sub('\\.[0-9]*$', '', rownames(combined_hsc))
#add tcgadf
TCGA_df
final_merged_counts <- na.omit(merge(TCGA_df, combined_hsc, by="row.names", all = TRUE))
final_merged_counts <- final_merged_counts %>% remove_rownames() %>% column_to_rownames(var = "Row.names")
final_merged_counts
write.table(final_merged_counts,file = "New_TCGA_HSC_merged_raw_counts_without_cas9.tsv",
            sep = "\t",quote = F)

#read coldata
col_data <- data.table::fread("../New_col_data_without_cas9.tsv")
col_data <- col_data %>% remove_rownames() %>% column_to_rownames(var = "ID")
col_data$Batch <- as.factor(col_data$Batch)


batch <- col_data$Batch
group <- ifelse(col_data$Type=="Control",0,1)


# generate pca plots 



##removing batch effects with RUVSeq and DESeq2
#run DESeq to determine control genes
dds <- DESeq2::DESeqDataSetFromMatrix (countData = as.data.frame(final_merged_counts), colData = col_data , design = ~ Type )
dds$Type <- relevel(dds$Type,ref = "Control")
#DE analysis before RUVSeq
dds <- DESeq(dds,test = "LRT",reduced = ~1)
res <- results(dds)

#run RUVSeq
library(RUVSeq)

set <- newSeqExpressionSet(counts(dds))
idx <- rowSums(counts(set) >5) >=2
set <- set[idx,]
set <- betweenLaneNormalization(set,which = "upper")
not.sig <- rownames(set)[which(res$pvalue>0.1)]
empirical <- rownames(set)[rownames(set) %in% not.sig]
#set number of factors
set <- RUVg(set,empirical, k=6)

#make plot of RUv factors
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:5) {
  stripchart(pData(set)[, i] ~ dds$SRA, vertical = TRUE, main = paste0("W", i))
  abline(h = 0)
}

#rerun DESEq with original dds 
#dds$W6 <- set_pdata$W_6
dds$W1 <- set$W_1
dds$W2 <- set$W_2
dds$W3 <- set$W_3
dds$W4 <- set$W_4
dds$W5 <- set$W_5
dds$W6 <- set$W_6
# dds$W7 <- set$W_7
design(dds) <- ~ W1 + W2 + W3 + W4+ W5+W6+ Type 
dds <- DESeq(dds)
dds <- dds[which(mcols(dds)$betaConv),]
resruv <- results(dds,alpha=0.05)
normalized_counts <- counts(dds,normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts)
result_after_ruv <- results(dds, contrast = c("Type", "Tumor", "Control")) 
result_after_ruv <- as.data.frame(result_after_ruv)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(result_after_ruv), mart= mart)
result_after_ruv$GeneName <- gene_IDs[match(rownames(result_after_ruv), gene_IDs[,1]),2]
result_after_ruv <- as.data.frame(result_after_ruv)
result_after_ruv$log2FoldChange <- as.numeric(result_after_ruv$log2FoldChange)
non_canonical_list <- read.table("~/Documents/AML/NonCanonical.txt", header = F)
raw_counts <- counts(dds,normalized=FALSE)

#plots
plot1 <- EnhancedVolcano::EnhancedVolcano(result_after_ruv[result_after_ruv$GeneName %in% non_canonical_list$V1,],
                                          lab = result_after_ruv[result_after_ruv$GeneName %in% non_canonical_list$V1,]$GeneName,
                                          x = 'log2FoldChange',
                                          y = 'padj',
                                          selectLab = c("GAPDH", "RPL18A","RPS10","RPL15","MRPL32"),
                                          title = 'NonCanonical RBPs',
                                          subtitle = '',
                                          caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.1'),
                                          pCutoff = 0.1,
                                          FCcutoff = 1,
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
                                          colConnectors = 'black')+
  ylim(c(0,7.5))+xlim(c(-5,5))
plot1
png("Batch_corrected/vol_plot_deseq.png",width = 35,height = 21,units = 'cm',res = 300)
plot(plot1)
dev.off()


##Remove batch effect from the counts matrix
vsd <- vst(dds,blind = FALSE)
mat <- assay(vsd)
mm <- model.matrix(~Type,colData(vsd))
mat2 <- limma::removeBatchEffect(mat,covariates = pData(set),design = mm)

##make PCA plot
p <- PCAtools::pca(mat2,metadata = col_data,removeVar = NULL)

#scree plot
scree_plot <- PCAtools::screeplot(p,axisLabSize = 18,titleLabSize = 22)
tiff("Batch_corrected/scree_plot_batch_corrected.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(scree_plot)
dev.off()
#biplot Plots 
pcat <- PCAtools::biplot(p,colby = 'Type',legendPosition = 'right',
                         lab=NULL)
tiff("Batch_corrected/PCA_tools_plot_batch_corrected_wo_labs.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(pcat)
dev.off()

#biplot with loading
pcat <- PCAtools::biplot(p, showLoadings = FALSE,
       labSize = 2, pointSize = 1.5, sizeLoadingsNames = 2,colby = "Type")
tiff("Batch_corrected/PCA_tools_plot_batch_corrected.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(pcat)
dev.off()

#replcae vsd with new matrix
assay(vsd) <- mat2
counts_batch_corrected <- assay(vsd)



#run wilcox test on batch corrected normalised counts 
conditions <- as.data.frame(t(col_data))[2,]
conditions <- factor(t(conditions))
count_norm <- mat2
pvalues <- sapply(1:nrow(count_norm), function(i){
  data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions)
  p <- wilcox.test(gene~conditions, data)$p.value
  return(p)})
fdr <- p.adjust(pvalues, method = "fdr")
# Calculate the fold-change for each gene
conditionsLevel <- levels(conditions)
dataCon1 <- count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2 <- count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
# Output results based on the FDR threshold 0.05
outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
rownames(outRst) <- rownames(count_norm)
outRst <- na.omit(outRst)
fdrThres <- 0.05
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(outRst), mart= mart)
outRst$GeneName <- gene_IDs[match(rownames(outRst), gene_IDs[,1]),2]


plot1 <- EnhancedVolcano::EnhancedVolcano(outRst[outRst$GeneName %in% non_canonical_list$V1,],
                                 lab = outRst[outRst$GeneName %in% non_canonical_list$V1,]$GeneName,
                                 x = 'log2foldChange',
                                 y = 'FDR',
                                 selectLab = c("GAPDH", "RPL18A","RPS10","RPL15","MRPL32"),
                                 title = 'NonCanonical RBPs',
                                 subtitle = '',
                                 caption = bquote(~Log[2]~ 'foldchange cutoff 0.5; p-value cutoff,0.05'),
                                 pCutoff = 0.05,
                                 FCcutoff = 0.5,
                                 pointSize = 1.0,
                                 labSize = 4.5,
                                 colAlpha = 1,
                                 legendPosition = 'right',
                                 legendLabSize = 8,
                                 legendIconSize = 2.0,
                                 drawConnectors = TRUE,
                                 widthConnectors = 0.95,
                                 labCol = 'black',
                                 labFace = 'bold',
                                 boxedLabels = TRUE,
                                 colConnectors = 'black')+
  xlim(c(-1,1))
tiff("Batch_corrected/vol_plot_wilcoxon_test.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(plot1)
dev.off()




#make violin plot
gene_id = "ENSG00000111640"
gene_name = "GAPDH"
GAPDH_all <- as.data.frame(t(normalized_counts[rownames(normalized_counts)=="ENSG00000111640",]))

GAPDH_all <- as.data.frame(mat2[rownames(mat2)=="ENSG00000111640",])
Control_samples <- rownames(col_data[col_data$Type=="Control",])
Tumor_samples <- rownames(col_data[col_data$Type!="Control",])
Control_samples = colnames(GAPDH_all[,which(colnames(mat2) %in% rownames(col_data[col_data$Type=="Control",]))])
Tumor_samples = colnames(GAPDH_all[,-which(colnames(mat2) %in% rownames(col_data[col_data$Type=="Control",]))])
# raw_GAPH <- t(final_merged_counts[rownames(final_merged_counts) == gene_id,])


GAPDH_all <- as.data.frame(GAPDH_all)
GAPDH_all$Type <- ifelse(rownames(GAPDH_all) %in% Control_samples, "Control", "Tumor")
#GAPDH_all$Raw <- raw_GAPH
colnames(GAPDH_all) <- c("Normalized", "Type")
GAPDH_all$Normalized[1:225] <- as.numeric(GAPDH_all$Normalized[1:225])

p1 <- ggplot(GAPDH_all[1:225,], aes(x = Type, y = Normalized, fill = Type))+
  geom_violin(trim = FALSE, width=0.4)+
  theme_cowplot()+
  geom_boxplot(width = 0.1, color = "grey", alpha=0.5)+
  theme(legend.position = "none",
        plot.title = element_text(size = 11))+
  ggtitle("GAPDH Normalized counts")+
  geom_jitter(height = 0, width = 0.1,alpha=0.7)+
  scale_fill_viridis_d()+
  scale_x_discrete(labels=c("Normal\n(n=74)","Tumor\n(n=151)"))+
  ylab("VST batch corrected counts") +
  xlab("Type")+
  geom_signif(comparisons = list(c("Control","Tumor")),
              map_signif_level = T,
              test = wilcox.test,
              y_position = 17.3, tip_length = 0.1, vjust = 0.2)+
  annotate("text",x = 0.58, xend = 1.9, y= 17, label = "Log2FC = 0.07")+
  annotate("text",x = 0.58, xend = 1.9, y= 16.8, label = "pvalue = 6.80e-26")
p1 
tiff("Batch_corrected/viol_plot.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(p1)
dev.off()

GAPDH_all %>% wilcox_effsize(Normalized ~ Type)
stat.test <- stat.test%>%add_xy_position(x="Type")

fold.change <- log(median(GAPDH_all[GAPDH_all$Type=="Tumor",1])/median(GAPDH_all[GAPDH_all$Type=="Control",1]),base=2)
fold.change

##----##
ggtitle("PCA with VST data")+
  theme_cowplot()




stat.test <- GAPDH_all %>%
  wilcox_test(Normalized ~ Type) %>%
  add_significance()
stat.test
fold.change <- log(mean(grslnd[grslnd$type=="Tumor",1])/mean(grslnd[grslnd$type=="Normal",1]),base=2)
fold.change <- log(median(grslnd[grslnd$type=="Tumor",1])/median(grslnd[grslnd$type=="Normal",1]),base=2)
#make PCA PLot
## create transformed values
vsd <- vst(dds,blind=FALSE)

pcaData <- plotPCA(vsd,intgroup=c("Type"),returnData=T)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
vst_pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Type,label=name)) +
  geom_point(size =1,alpha=0.4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")+
  theme_cowplot()+
  geom_text(hjust=1,size=1.5)
# GLM-PCA plots
tiff("Batch_corrected/plot_vst_with_label.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(vst_pca_plot)
dev.off()


#---------------------------------------------------------------------------------------------------#
#` SVA with DESeq2

dat <- counts(dds)















#---------------------------______---------------------------------____-----------------------------#
#tryin combatseq for batch correction
#batch correction with combatSeq
counts_without_batch <- sva::ComBat_seq(as.matrix(final_merged_counts),batch = batch,group = NULL,full_mod = TRUE)
#batch corrections with sva function
mod = model.matrix(~as.matrix(batch),)


`%nin%` = Negate(`%in%`)


