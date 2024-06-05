#--------------------------------------------------------------------------------------#
#GSEA for CLIPSeq common genes from union
#

#load the gene list
clip_gene_list <- data.table::fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Union/Gene_list_intersection.tsv",header = F)
clip_gene_list_int <- data.table::fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Intersection/Gene_list_all_intersection.tsv",header = F)
colnames(clip_gene_list) <- c("GeneID","GeneName")

venn.diagram(x=list(clip_gene_list$V1,
                    clip_gene_list_int$V1),
             category.names = c("Union","Int"),
             output=TRUE,
             filename = "venn.png", 
             col=c("#440154ff",'#fde725ff'),
             fill = c(alpha("#440154ff",0.3),  alpha('#fde725ff',0.3)),
             cex = 0.5,
             fontfamily = "sans",
             cat.cex = 0.3,
             cat.default.pos = "outer")

gapdh_gene_set <- data.frame("gs_name"=rep("GAPDH_TARGETS",length(unique(clip_gene_list$GeneID))),
                             "human_ensembl_gene"=unique(clip_gene_list$GeneID))
#add this to h_gene_set
gene_set <- rbind(h_gene_set,gapdh_gene_set)

#run gsea for high
generate_corr2 <- function(cts,cores=1,gene="ENSG00000111640") {
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
  return(corrMat)
}
corrMat <- generate_corr2(cts,cores = 6,gene = "ENSG00000111640")
corrMat$Spearman <- round(as.numeric(corrMat$Spearman),digits = 9)
corrMat <- corrMat[order(corrMat$Spearman,decreasing = TRUE),]
Spearman <- corrMat$Spearman
names(Spearman) <- corrMat$GeneID
pVals <- apply(as.matrix(Spearman), MARGIN = 1:2, n = length(corrDFFull[,1]), FUN = function(x, n) {
  stats::dt(abs(x)/sqrt((1-x^2)/(n-2)), df = 2)
})
pVals <- as.data.frame(pVals)


df <- as.data.frame(Spearman)
df$pval <- pVals[match(rownames(pVals),rownames(df)),1]
df$clip_target <- ifelse(rownames(df) %in% clip_gene_list$GeneID,
                         "Yes","No")
df$gene_name <- clip_gene_list[match(rownames(df),clip_gene_list$GeneID),]$GeneName

df_small <- df[df$clip_target=="Yes",]
df_small <- df_small[order(df_small$Spearman,decreasing = TRUE),]

gsea <- clusterProfiler::GSEA(Spearman,
                              minGSSize = 15,
                              maxGSSize = 500,
                              pAdjustMethod = "BH",
                              pvalueCutoff = 1,
                              seed = TRUE,
                              eps = 1e-100,
                              nPermSimple=1000,
                              by = "fgsea",
                              TERM2GENE = gene_set)
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
tiff("GSEA/Clip_targets/all_gapdh_gsea",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_hm)
dev.off()
erplot <- enrichplot::gseaplot2(gsea,geneSetID = gsea@result$ID[29],title = gsea@result$Description[29])
tiff("GSEA/Clip_targets/clip_targets_all_gseaplot.tiff",width = 35,height = 21,units = 'cm',res = 300)
erplot
dev.off()



heatmap <- ggplot(df_small[abs(df_small$Spearman)>0.3 & df_small$gene_name!="",],
                  aes(x=clip_target,y=reorder(gene_name,Spearman),
                      fill=Spearman))+
  geom_tile()+
  geom_text(aes(label=Spearman))+
  set_enrichplot_color(type = "fill",name="Correlation Factor")+
  theme_cowplot()+
  xlab("")+
  ylab("")+
  theme(axis.ticks.x =element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank())+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))
tiff("GSEA/Clip_targets/clip_correlations_high_gapdh.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(heatmap)
dev.off()


gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(df), mart= mart)
df$GeneName <- gene_IDs[match(rownames(df), gene_IDs[,1]),2]
write.table(df,file = "GSEA/Clip_targets/all_gapdh_correlations.tsv",quote = F,row.names = T,col.names = T,
            sep = "\t")
write.table(df_small,file = "GSEA/Clip_targets/all_gapdh_correlations_clip_targets.tsv",quote = F,row.names = T,col.names = T,
            sep = "\t")




