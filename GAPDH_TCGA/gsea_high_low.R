#create high vs low GAPDH expression
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(tidyverse)
#split gapdh expression from normalized counts used to make gene expression correlation
#tpm normalized
#old method, not to be used for gsea plots
#only use for correlation plots for individual genes
#replaced with using vst transformed counts from DESEQ2
high <- imp[imp$GAPDH>=median(imp$GAPDH),]
low <- imp[imp$GAPDH<median(imp$GAPDH),]
dim(high)
high_gapdh <- nm_tcga[,rownames(high)]
low_gapdh <- nm_tcga[,rownames(low)]
#make plot for high vs low
high_low <- imp
high_low$strata <- ifelse(imp$GAPDH>median(imp$GAPDH),"high","low") 

expr_plot <- ggplot(high_low, aes(x = strata, y = GAPDH, fill = strata))+
  geom_violin(trim = FALSE, width=0.4)+
  theme_cowplot()+
  geom_boxplot(width = 0.1, color = "grey", alpha=0.5)+
  theme(legend.position = "none",
        plot.title = element_text(size = 11))+
  geom_jitter(height = 0, width = 0.1,alpha=0.7)+
  scale_fill_viridis_d()+
  scale_x_discrete(labels=c("High\n(n=76)","Low\n(n=75)"))+
  ylab("Normalized Expression") +
  xlab("")+
  geom_signif(comparisons = list(c("high","low")),
              map_signif_level = T,
              test = wilcox.test,
              y_position = 5500, tip_length = 0.1, vjust = 0.2)
  # annotate("text",x = 0.58, xend = 1.9, y= 17, label = "Mean Log2FC = 0.1")+
  # annotate("text",x = 0.58, xend = 1.9, y= 16.8, label = "pvalue = 6.80e-26")
tiff("GSEA/high_vs_low.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(expr_plot)
dev.off()
  
GAPDH_all %>% wilcox_effsize(Normalized ~ Type)
stat.test <- stat.test%>%add_xy_position(x="Type")

## install.packages("devtools")
devtools::install_github("Bishop-Laboratory/correlationAnalyzeR")
library(correlationAnalyzeR)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(DOSE)
library(rstatix)


#deseq for tcga-laml data
nonZeroCount <- apply(TCGA_df_2, 1, function(row) {
  return(sum(row > 0))
})
keepInd <- which(nonZeroCount > (length(colnames(TCGA_df_2)) * .1))
TCGA_df <- TCGA_df_2[keepInd,]
dds <- DESeq2::DESeqDataSetFromMatrix(TCGA_df,colData = data.frame(sampleID=colnames(TCGA_df)),
                                      design = ~1)
inds <- rownames(TCGA_df)
geoMeansList <- parallel::mclapply(inds,FUN = function(ind){
  row <- TCGA_df[ind,]
  if (all(row == 0)) {
    0
  } else {
      exp(sum(log(row[row != 0]))/length(row))
  }
},mc.cores = 6)
geoMeansList <- unlist(geoMeansList)
dds <- DESeq2::estimateSizeFactors(dds,geoMeans=geoMeansList)
vsd <- DESeq2::vst(dds)
cts <- SummarizedExperiment::assay(vsd)
rownames(cts) <- sub('\\.[0-9]*$', '', rownames(cts))

#make high and low expression of gapdh data
imp <- cts["ENSG00000111640",]
imp <- as.data.frame(imp)
colnames(imp) <- c("GAPDH")
high <- rownames(imp)[which(imp$GAPDH>=median(imp$GAPDH))]
low <- rownames(imp)[which(imp$GAPDH<median(imp$GAPDH))]
dim(high)
high_gapdh <- cts[,high]
low_gapdh <- cts[,low]
#make plot for high vs low
high_low <- imp
high_low$strata <- ifelse(imp$GAPDH>median(imp$GAPDH),"high","low") 

attach(loadNamespace("enrichplot"),name = "enrichplot_all")


#make high vs low plot

expr_plot <- ggplot(high_low, aes(x = strata, y = GAPDH, fill = strata))+
  geom_violin(trim = FALSE, width=0.4)+
  theme_cowplot()+
  geom_boxplot(width = 0.1, color = "grey", alpha=0.5)+
  theme(legend.position = "none",
        plot.title = element_text(size = 11))+
  geom_jitter(height = 0, width = 0.1,alpha=0.7)+
  set_enrichplot_color(type = "color")+
  scale_x_discrete(labels=c("High\n(n=76)","Low\n(n=75)"))+
  ylab("Normalized Expression") +
  xlab("")+
  geom_signif(comparisons = list(c("high","low")),
              map_signif_level = T,
              test = wilcox.test,
              y_position = 18.2, tip_length = 0.1, vjust = 0.2)
# annotate("text",x = 0.58, xend = 1.9, y= 17, label = "Mean Log2FC = 0.1")+
# annotate("text",x = 0.58, xend = 1.9, y= 16.8, label = "pvalue = 6.80e-26")
tiff("GSEA/high_vs_low.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(expr_plot)
dev.off()




generate_corr <- function(cts,cores=1,gene="gene") {
  data <- as.data.frame(t(cts[-match(gene,rownames(cts)),]))
  gene <- as.data.frame(cts[match(gene,rownames(cts)),])
  # corrMat <- as.data.frame(t(cts)) %>% cor_mat("ENSG00000111640",method = "spearman")
  # corrPval <- corrMat %>% cor_get_pval()
  # corrMat <- rstatix::cor_mat(as.data.frame(t(cts)),y = as.matrix(gene),type =  "spearman")
  # corrMat <- cor.test(gene,data,method = "spearman")
  
  
  corrMat <- cor(gene,data,method = "spearman")
  corrMat <- data.frame(t(corrMat))
  corrMat$name <- rownames(corrMat)
  corrMat <- na.omit(corrMat)
  colnames(corrMat) <- c("Spearman","GeneID")
  # 
  #remove this condition to consider all correlation values
  #corrMat <- corrMat[abs(corrMat$Spearman)>0.3,]
  return(corrMat)
}

corrMat <- generate_corr(na.omit(high_gapdh),cores = 6,gene = "ENSG00000111640")
# corrMat$Spearman <- round(as.numeric(corrMat$Spearman),digits = 6)
corrMat <- corrMat[order(corrMat$Spearman,decreasing = TRUE),]
Spearman <- corrMat$Spearman
names(Spearman) <- corrMat$GeneID

#run GSEA analysis
h_gene_set <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name,human_ensembl_gene)
#run fgsea

gsea <- clusterProfiler::GSEA(Spearman,
                      minGSSize = 15,
                      maxGSSize = 500,
                      pAdjustMethod = "BH",
                      seed = TRUE,
                      eps = 1e-100,
                      nPermSimple=10000,
                      by = "fgsea",
                      TERM2GENE = h_gene_set)
gsea_neg <- clusterProfiler::GSEA(-1*Spearman[order(Spearman)],
                              minGSSize = 15,
                              maxGSSize = 500,
                              pAdjustMethod = "BH",
                              seed = TRUE,
                              eps = 1e-100,
                              nPermSimple=10000,
                              by = "fgsea",
                              TERM2GENE = h_gene_set)

gsea_result <- as.data.frame(gsea)
gsea_result_summary <- gsea_result[,c(1,3,5:8)]
gsea_result_summary$adjPvalue <- ifelse(gsea_result_summary$p.adjust <= 0.05, "significant", "non-significant")
cols <- c("non-significant" ="grey","significant"="red")
#general Plot
ggplot(gsea_result_summary, aes(reorder(ID, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")

#dotplot for upregulated
##make sorted dataframe
ggdata_up <- gsea_result_summary[gsea_result_summary$NES>0,]
ggdata_up$ID <- reorder(ggdata_up$ID,ggdata_up$NES,decreasing = F)


gsea_dotplot_up <- ggplot(ggdata_up[1:15,],aes(x=ID,y=NES,size=setSize,color=p.adjust,fill=p.adjust))+
  geom_point(shape=21)+
  scale_size(range = c(3,8))+
  scale_color_continuous(low='red',high='blue')+
  scale_fill_continuous(low='red',high='blue')+
  xlab('Gene set')+
  ylab('Enrichment Score')+
  labs(title = "High GAPDH (upregulated)")+
  theme_cowplot()+
  theme(
    legend.position = 'right',
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
    axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
    axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
    axis.title = element_text(size = 11, face = 'bold'),
    axis.title.y = element_text(size = 11, face = 'bold'),
    legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 12, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))+
  guides(size = guide_legend(order = 1))+
  coord_flip()
tiff("GSEA/high_gsea_up_dotplot.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_dotplot_up)
dev.off()

#dotplot for downregulated
##make sorted dataframe
ggdata_down <- gsea_result_summary[gsea_result_summary$NES<0,]
ggdata_down$ID <- reorder(ggdata_down$ID,ggdata_down$NES,decreasing = F)


gsea_dotplot_down <- ggplot(ggdata_down[1:15,],aes(x=ID,y=NES,size=setSize,color=p.adjust,fill=p.adjust))+
  geom_point(shape=21)+
  scale_size(range = c(3,8))+
  scale_color_continuous(low='red',high='blue')+
  scale_fill_continuous(low='red',high='blue')+
  xlab('Gene set')+
  ylab('Enrichment Score')+
  labs(title = "High GAPDH (downregulated)")+
  theme_cowplot()+
  theme(
    legend.position = 'right',
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
    axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
    axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
    axis.title = element_text(size = 11, face = 'bold'),
    axis.title.y = element_text(size = 11, face = 'bold'),
    legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 12, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))+
  guides(size = guide_legend(order = 1))+
  coord_flip()
tiff("GSEA/high_gsea_down_dotplot.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_dotplot_down)
dev.off()



glyco_gsea <- enrichplot::gseaplot2(gsea,geneSetID = gsea@result$ID[6],title = gsea@result$Description[6])
tiff("GSEA/high_glycolysis.tiff",width = 35,height = 21,units = 'cm',res = 300)
glyco_gsea
dev.off()

## Negative correlation genes
gsea_neg <- clusterProfiler::GSEA(-1*Spearman[order(Spearman)],
                                  minGSSize = 15,
                                  maxGSSize = 500,
                                  pAdjustMethod = "BH",
                                  seed = TRUE,
                                  eps = 1e-100,
                                  nPermSimple=10000,
                                  by = "fgsea",
                                  TERM2GENE = h_gene_set)

gsea_result_neg <- as.data.frame(gsea_neg)
gsea_result_summary_neg <- gsea_result_neg[,c(1,3,5:8)]
gsea_result_summary_neg$adjPvalue <- ifelse(gsea_result_summary_neg$p.adjust <= 0.05, "significant", "non-significant")
cols <- c("non-significant" ="grey","significant"="red")



#----------------------------------------------------------------------------------------##
# 
# # Low Gapdh Plots
# 
# corrMat <- generate_corr(low_gapdh,cores = 6,gene = "ENSG00000111640")
# corrMat$Spearman <- round(as.numeric(corrMat$Spearman),digits = 6)
# corrMat <- corrMat[order(corrMat$Spearman,decreasing = TRUE),]
# Spearman <- corrMat$Spearman
# names(Spearman) <- corrMat$GeneID
# 
# #run GSEA analysis
# #h_gene_set <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name,human_ensembl_gene)
# #run fgsea
# gsea <- clusterProfiler::GSEA(Spearman,
#                               minGSSize = 15,
#                               maxGSSize = 500,
#                               pAdjustMethod = "BH",
#                               seed = TRUE,
#                               eps = 1e-100,
#                               nPermSimple=1000,
#                               by = "fgsea",
#                               TERM2GENE = h_gene_set)
# gsea_result <- as.data.frame(gsea)
# gsea_result_summary <- gsea_result[,c(1,3,5:8)]
# gsea_result_summary$adjPvalue <- ifelse(gsea_result_summary$p.adjust <= 0.05, "significant", "non-significant")
# cols <- c("non-significant" ="grey","significant"="red")
# #general Plot
# ggplot(gsea_result_summary, aes(reorder(ID, NES), NES, fill = adjPvalue)) +
#   geom_col() +
#   scale_fill_manual(values = cols) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   coord_flip() +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title="Hallmark pathways Enrichment Score from GSEA")
# 
# #dotplot for upregulated
# ##make sorted dataframe
# ggdata_up <- gsea_result_summary[gsea_result_summary$NES>0,]
# ggdata_up$ID <- reorder(ggdata_up$ID,ggdata_up$NES,decreasing = F)
# 
# 
# gsea_dotplot_up <- ggplot(ggdata_up[1:15,],aes(x=ID,y=NES,size=setSize,color=p.adjust,fill=p.adjust))+
#   geom_point(shape=21)+
#   scale_size(range = c(3,8))+
#   scale_color_continuous(low='red',high='blue')+
#   scale_fill_continuous(low='red',high='blue')+
#   xlab('Gene set')+
#   ylab('Enrichment Score')+
#   labs(title = "Low GAPDH (upregulated)")+
#   theme_cowplot()+
#   theme(
#     legend.position = 'right',
#     plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
#     axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
#     axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
#     axis.title = element_text(size = 11, face = 'bold'),
#     axis.title.y = element_text(size = 11, face = 'bold'),
#     legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
#     legend.text = element_text(size = 12, face = "bold"), # Text size
#     title = element_text(size = 14, face = "bold"))+
#   guides(size = guide_legend(order = 1))+
#   coord_flip()
# tiff("GSEA/low_gsea_up_dotplot.tiff",width = 35,height = 21,units = 'cm',res = 300)
# plot(gsea_dotplot_up)
# dev.off()
# 
# #dotplot for downregulated
# ##make sorted dataframe
# ggdata_down <- gsea_result_summary[gsea_result_summary$NES<0,]
# ggdata_down$ID <- reorder(ggdata_down$ID,ggdata_down$NES,decreasing = F)
# 
# 
# gsea_dotplot_down <- ggplot(ggdata_down[1:15,],aes(x=ID,y=NES,size=setSize,color=p.adjust,fill=p.adjust))+
#   geom_point(shape=21)+
#   scale_size(range = c(3,8))+
#   scale_color_continuous(low='red',high='blue')+
#   scale_fill_continuous(low='red',high='blue')+
#   xlab('Gene set')+
#   ylab('Enrichment Score')+
#   labs(title = "Low GAPDH (downregulated)")+
#   theme_cowplot()+
#   theme(
#     legend.position = 'right',
#     plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
#     axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
#     axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
#     axis.title = element_text(size = 11, face = 'bold'),
#     axis.title.y = element_text(size = 11, face = 'bold'),
#     legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
#     legend.text = element_text(size = 12, face = "bold"), # Text size
#     title = element_text(size = 14, face = "bold"))+
#   guides(size = guide_legend(order = 1))+
#   coord_flip()
# tiff("GSEA/low_gsea_down_dotplot.tiff",width = 35,height = 21,units = 'cm',res = 300)
# plot(gsea_dotplot_down)
# dev.off()
# 
# 
# 
# glyco_gsea <- enrichplot::gseaplot2(gsea,geneSetID = gsea@result$ID[6],title = gsea@result$Description[6])
# tiff("GSEA/low_glycolysis.tiff",width = 35,height = 21,units = 'cm',res = 300)
# glyco_gsea
# dev.off()
# 
# # low gapdh, high expresssion genes
# gsea_neg <- clusterProfiler::GSEA(-1*Spearman[order(Spearman)],
#                                   minGSSize = 15,
#                                   maxGSSize = 500,
#                                   pAdjustMethod = "BH",
#                                   seed = TRUE,
#                                   eps = 1e-100,
#                                   nPermSimple=10000,
#                                   by = "fgsea",
#                                   TERM2GENE = h_gene_set)
# 
# gsea_result <- as.data.frame(gsea_neg)
# gsea_result_summary <- gsea_result[,c(1,3,5:8)]
# gsea_result_summary$adjPvalue <- ifelse(gsea_result_summary$p.adjust <= 0.05, "significant", "non-significant")
# cols <- c("non-significant" ="grey","significant"="red")
# 
# 
# #dotplot for upregulated
# ##make sorted dataframe
# ggdata_up <- gsea_result_summary[gsea_result_summary$NES>0,]
# ggdata_up$ID <- reorder(ggdata_up$ID,ggdata_up$NES,decreasing = F)
# 
# 
# gsea_dotplot_up <- ggplot(ggdata_up[1:15,],aes(x=ID,y=NES,size=setSize,color=p.adjust,fill=p.adjust))+
#   geom_point(shape=21)+
#   scale_size(range = c(3,8))+
#   scale_color_continuous(low='red',high='blue')+
#   scale_fill_continuous(low='red',high='blue')+
#   xlab('Gene set')+
#   ylab('Enrichment Score')+
#   labs(title = "Low GAPDH (upregulated)")+
#   theme_cowplot()+
#   theme(
#     legend.position = 'right',
#     plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
#     axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
#     axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
#     axis.title = element_text(size = 11, face = 'bold'),
#     axis.title.y = element_text(size = 11, face = 'bold'),
#     legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
#     legend.text = element_text(size = 12, face = "bold"), # Text size
#     title = element_text(size = 14, face = "bold"))+
#   guides(size = guide_legend(order = 1))+
#   coord_flip()
# tiff("GSEA/low_gsea_up_dotplot.tiff",width = 35,height = 21,units = 'cm',res = 300)
# plot(gsea_dotplot_up)
# dev.off()
# 
# #dotplot for downregulated
# ##make sorted dataframe
# ggdata_down <- gsea_result_summary[gsea_result_summary$NES<0,]
# ggdata_down$ID <- reorder(ggdata_down$ID,ggdata_down$NES,decreasing = F)
# 
# 
# gsea_dotplot_down <- ggplot(ggdata_down[1:15,],aes(x=ID,y=NES,size=setSize,color=p.adjust,fill=p.adjust))+
#   geom_point(shape=21)+
#   scale_size(range = c(3,8))+
#   scale_color_continuous(low='red',high='blue')+
#   scale_fill_continuous(low='red',high='blue')+
#   xlab('Gene set')+
#   ylab('Enrichment Score')+
#   labs(title = "Low GAPDH (downregulated)")+
#   theme_cowplot()+
#   theme(
#     legend.position = 'right',
#     plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
#     axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
#     axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
#     axis.title = element_text(size = 11, face = 'bold'),
#     axis.title.y = element_text(size = 11, face = 'bold'),
#     legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
#     legend.text = element_text(size = 12, face = "bold"), # Text size
#     title = element_text(size = 14, face = "bold"))+
#   guides(size = guide_legend(order = 1))+
#   coord_flip()
# tiff("GSEA/low_gsea_down_dotplot.tiff",width = 35,height = 21,units = 'cm',res = 300)
# plot(gsea_dotplot_down)
# dev.off()
# 
# 
# 

#---------------------------------------------------------------------------------_____###
# All samples(high + low)

corrMat <- generate_corr(cts,cores = 6,gene = "ENSG00000111640")
corrMat$Spearman <- round(as.numeric(corrMat$Spearman),digits = 6)
corrMat <- corrMat[order(corrMat$Spearman,decreasing = TRUE),]
Spearman <- corrMat$Spearman
names(Spearman) <- corrMat$GeneID

#run GSEA analysis
##Load hallmark data
#h_gene_set <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name,human_ensembl_gene)
#run fgsea
gsea <- clusterProfiler::GSEA(Spearman,
                              minGSSize = 15,
                              maxGSSize = 500,
                              pAdjustMethod = "BH",
                              seed = TRUE,
                              eps = 1e-100,
                              nPermSimple=1000,
                              by = "fgsea",
                              TERM2GENE = h_gene_set)
gsea_result <- as.data.frame(gsea)
gsea_result_summary <- gsea_result[,c(1,3,5:8)]
gsea_result_summary$adjPvalue <- ifelse(gsea_result_summary$p.adjust <= 0.05, "significant", "non-significant")
cols <- c("non-significant" ="grey","significant"="red")

#general Plot
gsea_hm <- ggplot(gsea_result_summary, aes(reorder(ID, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")
tiff("GSEA/all_pathway.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_hm)
dev.off()

#dotplot for upregulated
##make sorted dataframe
ggdata_up <- gsea_result_summary[gsea_result_summary$NES>0,]
ggdata_up$ID <- reorder(ggdata_up$ID,ggdata_up$NES,decreasing = F)


gsea_dotplot_up <- ggplot(ggdata_up[1:15,],aes(x=ID,y=NES,size=setSize,color=p.adjust,fill=p.adjust))+
  geom_point(shape=21)+
  scale_size(range = c(3,8))+
  scale_color_continuous(low='red',high='blue')+
  scale_fill_continuous(low='red',high='blue')+
  xlab('Gene set')+
  ylab('Enrichment Score')+
  labs(title = "Summary dotplot (upregulated)")+
  theme_cowplot()+
  theme(
    legend.position = 'right',
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
    axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
    axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
    axis.title = element_text(size = 11, face = 'bold'),
    axis.title.y = element_text(size = 11, face = 'bold'),
    legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 12, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))+
  guides(size = guide_legend(order = 1))+
  coord_flip()
tiff("GSEA/all_gsea_up_dotplot.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_dotplot_up)
dev.off()

#dotplot for downregulated
##make sorted dataframe
ggdata_down <- gsea_result_summary[gsea_result_summary$NES<0,]
ggdata_down$ID <- reorder(ggdata_down$ID,ggdata_down$NES,decreasing = F)


gsea_dotplot_down <- ggplot(ggdata_down[1:15,],aes(x=ID,y=NES,size=setSize,color=p.adjust,fill=p.adjust))+
  geom_point(shape=21)+
  scale_size(range = c(3,8))+
  scale_color_continuous(low='red',high='blue')+
  scale_fill_continuous(low='red',high='blue')+
  xlab('Gene set')+
  ylab('Enrichment Score')+
  labs(title = "Summary dotplot (downregulated)")+
  theme_cowplot()+
  theme(
    legend.position = 'right',
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1, hjust = 0.5, color = "black"),
    axis.text.x = element_text(angle = 0, size = 11, face = 'bold', hjust = 0.5, color = "black"),
    axis.text.y = element_text(angle = 0, size = 11, face = 'bold', vjust = 0.5, color = "black"),
    axis.title = element_text(size = 11, face = 'bold'),
    axis.title.y = element_text(size = 11, face = 'bold'),
    legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 12, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold"))+
  guides(size = guide_legend(order = 1))+
  coord_flip()
tiff("GSEA/all_gsea_down_dotplot.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_dotplot_down)
dev.off()



glyco_gsea <- enrichplot::gseaplot2(gsea,geneSetID = gsea@result$ID[5],title = gsea@result$Description[5])
tiff("GSEA/all_glycolysis.tiff",width = 35,height = 21,units = 'cm',res = 300)
glyco_gsea
dev.off()




