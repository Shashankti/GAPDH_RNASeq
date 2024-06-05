library(tidyverse)
library(ggplot2)
library(ggrepel)
library(DEGreport)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(biomaRt)
library(viridis)
library(genefilter)
library(ComplexHeatmap)
library(data.table)
library(cowplot)


#Make lfc shrinkage for shNorm vs ctNorm
#res_shNorm_vs_ctNorm <- lfcShrink(dds_2f,coef = "type_shRNA_vs_parental",type = "apeglm")
res_shNorm_vs_ctNorm <- lfcShrink(dds_2f,contrast = norm_shRNA - norm_parental,type = "ashr")
GeneList <- 


# Generate headmap for Norm ct vs norm shrna
vsd <- vst(dds_2f,blind=FALSE)
rld <- rlog(dds_2f,blind=F)
head(assay(vsd),3)
select <- order(rowMeans(counts(dds_2f,normalized=TRUE)),
                decreasing=TRUE)[1:40]
df <- as.data.frame(colData(dds_2f)[,c("condition","type")])

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df,scale = "row")
#$PCA Plots

pcaData <- plotPCA(vsd,intgroup=c("type","condition"),returnData=T)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
vst_pca_plot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = type, shape = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA plot of RNA-Seq Samples")+
  theme_cowplot()
#save the plot
ggsave("Plots/PCA_Plot.jpeg",vst_pca_plot,width = 20, height = 10, units = 'in',dpi = 300)
plot(vst_pca_plot)
#save the output plot as jpeg
jpeg("Plots/Samples_PCA_Plot.jpeg",width = 35,height = 21,units = 'cm',res = 300)
grid.arrange(vst_pca_plot)
dev.off()
# GLM-PCA plots
library("glmpca")
gpca <- glmpca(counts(dds_2f), L=2)
gpca.dat <- gpca$factors
gpca.dat$type <- dds_2f$type
gpca.dat$condition <- dds_2f$condition

gpca_plot <- ggplot(gpca.dat, aes(x = dim1, y = dim2, color = type, shape = condition)) +
  geom_point(size =3) + coord_fixed() + ggtitle("Generalized PCA plot of RNA-Seq Samples")+
  theme_cowplot()
jpeg("Plots/Samples_GPCA_Plot.jpeg",width = 35,height = 21,units = 'cm',res = 300)
grid.arrange(gpca_plot)
dev.off()


# Define Samptype  info in META df
sampletype <- factor(c(rep("parental",6),rep("shRNA",6)))
meta <- data.frame(sampletype,row.names = rownames(colData(dds_2f)))


heat.colors <- brewer.pal(12, "Set3")


# Measure and plot Sample distances
sampleDists <- dist(t(assay(rld)))
as.matrix(sampleDists)[1:2,1:2]
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition,rld$type,sep = "-")
colnames(sampleDistMatrix) <- NULL
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
gplots::heatmap.2( sampleDistMatrix, trace="none", col=colours)


#Geneclustering with rld object
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
#Better looking heatmaps
gplots::heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[
             colData(rld)$type ] )



# Heat map with VSD data
vsd_mat <- assay(vsd)
vsd_mat <- as.data.frame(vsd_mat)
vsd_mat$GeneID <- rownames(vsd_mat)

sigGenes <- res_shNorm_vs_ctNorm_2[which(res_shNorm_vs_ctNorm_2$padj <= 0.1 & abs(res_shNorm_vs_ctNorm_2$log2FoldChange)>0.2),]$GeneID
vsd_mat <- vsd_mat[vsd_mat$GeneID %in% sigGenes,]
normalized_counts_mat <- normalized_counts[normalized_counts$GeneID %in% sigGenes,]
#shink
vsd_wide <- vsd_mat
vsd_long <- reshape2::melt(vsd_mat,id.vars=c("GeneID"))
norm_long <- reshape2::melt(normalized_counts_mat,vars=c("GeneID"))
vsd_long$FillScaled <- as.numeric(scale(vsd_long$value))
norm_long$FillScaled <- as.numeric(scale(norm_long$value))

heatmap <- ggplot(vsd_long,aes(x=variable,y=GeneID,fill=FillScaled))+ geom_raster()+scale_fill_viridis()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1),axis.text.y = element_text(angle = 45,hjust = 1),axis.ticks.y = element_blank())


heatmap1 <- ggplot(norm_long,aes(x=variable,y=GeneID,fill=FillScaled))+ geom_raster()+scale_fill_viridis()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1),axis.text.y = element_text(angle = 45,hjust = 1),axis.ticks.y = element_blank())+
  scale_color_gradient2(midpoint = 0)
#normalised counts heatmap
pheatmap(normalized_counts_mat,
         color = heat.colors,
         show_rownames = T,
         annotation = meta,
         scale = "row")

#Volcano Plot
ggplot()
DEGreport::degVolcano(
  data.frame(res_shNorm_vs_ctNorm[,c("log2FoldChange","padj")]))
DEGreport::degVolcano(Gal_df[,c("log2FoldChange","qvalue")])







# Log Fold change shNorm vs ctNorm
N <- 20
top_genes <- res_shNorm_vs_ctNorm_2$GeneID[1:N]
res_shNorm_vs_ctNorm_2 %>%
  mutate(Label = ifelse(GeneID %in% top_genes,hgnc_symbol,"")) %>%
  ggplot(aes(x= log2FoldChange,y=-log10(padj),label=Label))+ geom_point(alpha=0.4 + geom_text(col="blue"))

ggplot(res_shNorm_vs_ctNorm_2[1:20,]) + geom_point(aes(x=log2FoldChange,y= -log10(padj),label=hgnc_symbol),alpha=0.4)


# Plot for one gene
d <- plotCounts(dds_2f,gene =which.min(res_shNorm_vs_ctNorm_2$padj) ,intgroup = "condition",returnData = TRUE)
ggplot(d, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("MOV10") +
  theme(plot.title = element_text(hjust = 0.5))







#Get norm counts again
normalized_counts <- normalized_counts[which(isUnique(normalized_counts$GeneID)),]
rownames(normalized_counts) <- normalized_counts$GeneID


##Set cutoffs
padj.cutoff <- 0.1
lfc.cutoff <- 0.5

threshold <- res_shNorm_vs_ctNorm$padj < padj.cutoff & abs(res_shNorm_vs_ctNorm$log2FoldChange) > lfc.cutoff
length(which(threshold))
res_shNorm_vs_ctNorm$threshold <- threshold
sigOE <- data.frame(subset(res_shNorm_vs_ctNorm,threshold==TRUE))
rownames(normalized_counts) <- normalized_counts$GeneID
sigOE_ordered <- sigOE[order(sigOE$padj),]
sigOE_ordered <- normalized_counts[rownames(sigOE_ordered),]
melted_sigOE <- data.frame(reshape::melt(sigOE_ordered))
melted_sigOE <- na.omit(melted_sigOE)
meta$samplename <- rownames(meta)
melted_sigOE$type <- meta[match(melted_sigOE$variable,meta$samplename),1]
#colnames(melted_sigOE) <-c("GeneID","condition","value","type")

ggplot(melted_sigOE)+
  geom_point(aes(x=GeneName,y=value,color=type))+
  scale_y_log10() + 
  xlab("Genes")+
  ylab("Normalized Counts")+
  ggtitle("Significant DE genes")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))

for(i in 1:length(melted_sigOE$GeneID)){
  melted_sigOE$FillScaled[i] <- (melted_sigOE$value[i] - mean(melted_sigOE$value))/sd(melted_sigOE$value) 
}

ggplot(melted_sigOE,aes(x=variable,y=GeneName,fill=FillScaled))+ 
  geom_raster()+
  scale_fill_viridis_c(name="Normalized Counts")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1),axis.text.y = element_text(angle = 45,hjust = 1,size = 7),axis.ticks.y = element_blank())
sigOE_ordered <- sigOE_ordered %>% mutate_at(c(2:13),as.numeric)
sig_mat <- as.matrix(sigOE_ordered[,c(2:13)])
`class<-` (sig_mat,'numeric')
meta$samplename <- NULL
pheatmap(na.omit(sig_mat),
         show_rownames = T,
         annotation_col = meta,
         color = brewer.pal(3,"RdBu"),
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         main = "Normalized counts for Signigicant DE genes")

sig_mat_scaled <-t(scale(t(sigOE_ordered[,c(2:13)])))
rownames(sig_mat_scaled) <- sigOE_ordered$GeneName

ht_all <-ComplexHeatmap::Heatmap(sig_mat_scaled,
                                  column_names_side = "bottom",
                                  cluster_rows = T,
                                  cluster_columns = F,
                                  column_names_rot = 45,
                                  column_names_gp = gpar(fontsize=7),
                                  row_names_gp = gpar(fontsize=7),
                                  row_names_rot = -15,
                                 name = "z-score")

# up-regulated genes
padj.cutoff <- 0.1
lfc.cutoff <- 0.5

threshold <- res_shNorm_vs_ctNorm$padj < padj.cutoff & res_shNorm_vs_ctNorm$log2FoldChange > lfc.cutoff
length(which(threshold))
res_shNorm_vs_ctNorm$threshold <- threshold
sigOE_up <- data.frame(subset(res_shNorm_vs_ctNorm,threshold==TRUE))

sigOE_ordered_up <- sigOE_up[order(sigOE_up$padj),]
sigOE_ordered_up <- normalized_counts[rownames(sigOE_ordered_up),]
melted_sigOE_up <- data.frame(reshape::melt(sigOE_ordered_up))
meta$samplename <- rownames(meta)
melted_sigOE_up$type <- meta[match(melted_sigOE_up$variable,meta$samplename),1]

ggplot(melted_sigOE_up)+
  geom_point(aes(x=GeneName,y=value,color=type))+
  scale_y_log10() + 
  xlab("Genes")+
  ylab("Normalized Counts")+
  ggtitle("Up regulated Significant DE genes")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))

for(i in 1:length(melted_sigOE_up$GeneID)){
  melted_sigOE_up$FillScaled[i] <- (melted_sigOE_up$value[i] - mean(melted_sigOE_up$value))/sd(melted_sigOE_up$value) 
}

ggplot(melted_sigOE_up,aes(x=variable,y=GeneName,fill=FillScaled))+ 
  geom_raster()+
  scale_fill_viridis_c(name="Normalized Counts")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1),axis.text.y = element_text(angle = 45,hjust = 1,size = 7),axis.ticks.y = element_blank())
sigOE_ordered_up <- sigOE_ordered_up %>% mutate_at(c(2:13),as.numeric)
sig_mat_up <- as.matrix(sigOE_ordered_up[,c(2:13)])
`class<-` (sig_mat_up,'numeric')
rownames(sig_mat_up) <- 
sigOE_ordered_up$GeneName
meta$samplename <- NULL
pheatmap(sig_mat_up[,c(1,2,7,10,3,4,8,11,5,6,9,12)],
         show_rownames = T,
         annotation_col = meta,
         color = brewer.pal(3,"RdBu"),
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         main = "Up-regulated Significant DE genes")

##################3
#Down regulated genes

padj.cutoff <- 0.1
lfc.cutoff <- -0.5

threshold <- res_shNorm_vs_ctNorm$padj < padj.cutoff & res_shNorm_vs_ctNorm$log2FoldChange < lfc.cutoff
length(which(threshold))
res_shNorm_vs_ctNorm$threshold <- threshold
sigOE_down <- data.frame(subset(res_shNorm_vs_ctNorm,threshold==TRUE))

threshold <-deRes$qvalue < padj.cutoff & deRes$log2FoldChange < lfc.cutoff
length(which(threshold))
deRes$threshold <- threshold
sigOE_down <- data.frame(subset(deRes,threshold==TRUE))

sigOE_ordered_down <- sigOE_down[order(sigOE_down$padj),]
sigOE_ordered_down <- normalized_counts[rownames(sigOE_ordered_down),]
melted_sigOE_down <- data.frame(reshape::melt(sigOE_ordered_down))
meta$samplename <- rownames(meta)
melted_sigOE_down$type <- meta[match(melted_sigOE_down$variable,meta$samplename),1]

ggplot(melted_sigOE_down)+
  geom_point(aes(x=GeneName,y=value,color=type))+
  scale_y_log10() + 
  xlab("Genes")+
  ylab("Normalized Counts")+
  ggtitle("Down regulatedSignificant DE genes")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))

for(i in 1:length(melted_sigOE_down$GeneID)){
  melted_sigOE_down$FillScaled[i] <- (melted_sigOE_down$value[i] - mean(melted_sigOE_down$value))/sd(melted_sigOE_down$value) 
}

ggplot(melted_sigOE_down,aes(x=variable,y=GeneName,fill=FillScaled))+ 
  geom_raster()+
  scale_fill_viridis_c(name="Normalized Counts")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1),axis.text.y = element_text(angle = 45,hjust = 1,size = 7),axis.ticks.y = element_blank())
sigOE_ordered_down <- sigOE_ordered_down %>% mutate_at(c(2:13),as.numeric)
sig_mat_down <- as.matrix(sigOE_ordered_down[,c(2:13)])
`class<-` (sig_mat_down,'numeric')
rownames(sig_mat_down) <- 
  sigOE_ordered_down$GeneName

meta$samplename <- NULL
pheatmap(sig_mat_down[,c(1,2,7,10,3,4,8,11,5,6,9,12)],
         show_rownames = F,
         annotation_col = meta,
         color = brewer.pal(3,"RdBu"),
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         main = "Down regulated Significant DE genes")



pheatmap(res_shNorm_vs_ctNorm$log2FoldChange,
         annotation_col = meta,
         color = brewer.pal(3,"Blues"))












norm_OEsig <- normalized_counts[rownames(sigOE),]
norm_OEsig$GeneID <- NULL
norm_OEsig$GeneName <- NULL
norm_OEsig <- as.matrix(norm_OEsig)
pheatmap(norm_OEsig, color = heat.colors,
         breaks = )
pheatmap(norm_OEsig, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE,scale = "row")

hm.expr <- plot_heatmap.expr(data = normalized_counts, 
                             annot = colData_2f,
                             sample.id = "sample_id",
                             annot.flags = c("condition", "type"),
                             ct.table.id.type = "ENSEMBL",
                             gene.id.type = "SYMBOL",
                             gene.names = NULL,
                             gene.count = 10,
                             title = "RVA Heatmap",
                             fill = "CPM",
                             baseline.flag = "control",
                             baseline.val = "0",
                             plot.save.to = NULL,
                             input.type = "count")
colData_2f <- as.data.frame(colData_2f) %>% rownames_to_column("sample_id")



## Trying reporting tools library to generate html report
dat  <- counts(dds_2f, normalized = TRUE)
idx  <- rowMeans(dat) > 1
dat  <- dat[idx, ]
mod  <- model.matrix(~ type, colData(dds_2f))
mod0 <- model.matrix(~   condition, colData(dds_2f))
svseq <- svaseq(dat, mod, mod0, n.sv = 2)

par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[, i] ~ dds_2f$condition, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
}

ddssva <- dds_2f
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + condition + type + type:condition
ddssva <- DESeq(ddssva)
res_sva_norm <- results(ddssva, contrast = norm_shRNA - norm_parental)
