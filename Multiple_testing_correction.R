#! For IHW and alternate methods to get better padj values

library(IHW)
library(dplyr)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)





#Load the Deseq results
deRes <- as.data.frame(res_shNorm_vs_ctNorm_Unshrunk)
#Plot of the pvalue distribution
ggplot(as.data.frame(res_pLow_vs_pNorm_2), aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0,fill="azure4")+
  theme_cowplot()+
  ggtitle("p-value distribution without correction")
#Info regarding metadata
metadata(res_shNorm_vs_ctNorm_Unshrunk)$filterThreshold

#Plot for independent filtering
plot(metadata(res_shLow_vs_ctLow_unshrunk)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res_shLow_vs_ctLow_unshrunk)$lo.fit, col="red")
abline(v=metadata(res_shLow_vs_ctLow_unshrunk)$filterTheta)

# Check the filtering criteria

plot(res_shGal_vs_ctGal_unshrunk$baseMean+1, -log10(res_shGal_vs_ctGal_unshrunk$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
use <- res_shGal_vs_ctGal_unshrunk$baseMean > metadata(res_shGal_vs_ctGal_unshrunk)$filterThreshold
h1 <- hist(res_shGal_vs_ctGal_unshrunk$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res_shGal_vs_ctGal_unshrunk$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

#Apply fdrtool to value
deRes.fixed <- fdrtool::fdrtool(deRes$stat,statistic = "normal")
deRes$qvalue <- p.adjust(deRes.fixed$pval,method = "BH")
deRes <- arrange(deRes,qvalue) %>% rownames_to_column("GeneID")

p.adj <- p.adjust(res_shNorm_vs_ctNorm_Unshrunk$pvalue,method = "BH")


#save output as csv

write.csv(deRes,file = "Analysis/res_shNorm_vs_ctNorm")
# Making heatmaps with different adjusted p-values

padj.cutoff <- 0.05
lfc.cutoff <- 0.5

threshold <- deRes$qvalue < padj.cutoff & abs(deRes$log2FoldChange) > lfc.cutoff
length(which(threshold))
deRes$threshold <- threshold
sigOE_qValue <- data.frame(subset(deRes,threshold==TRUE))


sigOE_ordered_qValue <- sigOE_qValue[order(sigOE_qValue$qvalue),]
sigOE_ordered_qValue <- normalized_counts[rownames(sigOE_ordered_qValue),]
sigOE_ordered_qValue$qvalue <- sigOE_qValue[match(sigOE_ordered_qValue$GeneID,rownames(sigOE_qValue)),7]
melted_sigOE_qValue <- data.frame(reshape::melt(sigOE_ordered_qValue[,c(1:14)]))
meta$samplename <- rownames(meta)

melted_sigOE_qValue$type <- meta[match(melted_sigOE_qValue$variable,meta$samplename),1]
melted_sigOE_qValue$lfc <- sigOE_qValue[match(melted_sigOE_qValue$GeneID,rownames(sigOE_qValue)),2]
melted_sigOE_qValue$qvalue <- sigOE_qValue[match(melted_sigOE_qValue$GeneID,rownames(sigOE_qValue)),7]
unique(melted_sigOE_qValue$qvalue)[1:20]

# up-regulated genes
ggplot(melted_sigOE_qValue[melted_sigOE_qValue$qvalue %in% unique(melted_sigOE_qValue$qvalue)[1:30] & melted_sigOE_qValue$lfc<0,])+
  geom_point(aes(x=GeneName,y=value,color=type))+
  scale_y_log10() + 
  xlab("Genes")+
  ylab("Normalized Counts")+
  ggtitle("Top 20 Down-regulated Significant DE genes")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))

for(i in 1:length(melted_sigOE_qValue$GeneID)){
  melted_sigOE_qValue$FillScaled[i] <- (melted_sigOE_qValue$value[i] - mean(melted_sigOE_qValue$value))/sd(melted_sigOE_qValue$value) 
}

ggplot(melted_sigOE_qValue,aes(x=variable,y=GeneName,fill=FillScaled))+ 
  geom_raster()+
  scale_fill_viridis_c(name="Normalized Counts")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1),axis.text.y = element_text(angle = 45,hjust = 1,size = 7),axis.ticks.y = element_blank())
sigOE_ordered_qValue <- sigOE_ordered_qValue %>% mutate_at(c(2:13),as.numeric)
#scale the df
mat_scaled <- t(scale(t(sigOE_ordered_qValue[,c(2:13)])))
rownames(mat_scaled) <- sigOE_ordered_qValue$GeneName

sig_mat_qValue <- as.matrix(sigOE_ordered_qValue[,c(2:13)])
`class<-` (sig_mat_qValue,'numeric')
rownames(sig_mat_qValue) <- 
  sigOE_ordered_qValue$GeneName
meta$samplename <- NULL
pheatmap(sig_mat_qValue,
         show_rownames = T,
         annotation_col = meta,
         scale = "row",
         cluster_cols = F,
         cluster_rows = T,
         main = "Significant DE genes")
#2
gplots::heatmap.2( sig_mat_qValue, scale="row", 
                   trace="none", dendrogram="column", 
                   col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
                   ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[
                     colData(dds_2f)$type ] )
#2
ht <- ComplexHeatmap::Heatmap(mat_scaled,
                       column_names_side = "bottom",
                        cluster_rows = T,
                        cluster_columns = F,
                        column_names_rot = 45,
                       column_names_gp = gpar(fontsize=7),
                       row_names_gp = gpar(fontsize=7),
                       row_names_rot = -15)



#output csv file
sigOE_qValue$GeneName <- sigOE_ordered_qValue[match(rownames(sigOE_qValue),sigOE_ordered_qValue$GeneID),14]
write.csv(sigOE_qValue,"Analysis/res_shNorm_vs_ctNorm_sig.csv")

#Making plot with melted Df.

only_norm <- sigOE_qValue[order(sigOE_qValue$qvalue),]
only_norm <- normalized_counts[rownames(only_norm),c(1,2,3,8,11,14)]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
only_norm_mat <- as.matrix(only_norm[,c(2,3,4,5)]) 
only_norm_mat_scaled <- t(scale(t(only_norm[,c(2:5)])))
rownames(only_norm_mat_scaled) <- only_norm$GeneName
ht_norm <-ComplexHeatmap::Heatmap(only_norm_mat_scaled,
                        cluster_rows = T,
                        cluster_columns = F,
                        column_names_rot = 45,
                        column_names_gp = gpar(fontsize=7),
                        row_names_gp = gpar(fontsize=4),
                        row_names_rot = 0,
                        name = "Z-score",width = 4*unit(25,'mm'),
                        column_title = "controlNorm_vs_shNorm")
ht_norm
tiff("Plots/shRNA_vs_ctNorm_heatmap.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_norm)
dev.off()



p <- plotly::plot_ly(x=colnames(only_norm_mat_scaled),y=rownames(only_norm_mat_scaled),z = only_norm_mat_scaled,type = "heatmap")

melt_norm <- data.frame(reshape::melt(only_norm_mat_scaled))
p1 <- ggplot(melt_norm,aes(X2,X1),col="RdYlBu")+
  geom_raster(aes(fill=value))+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 65,hjust = 1),axis.text.y = element_text(angle = 45,hjust = 1,size = 7),axis.ticks.y = element_blank())

  
plotly::ggplotly(p1)


# making volcano plot
library(EnhancedVolcano)

deRes$GeneName <- normalized_counts[match(rownames(deRes),normalized_counts$GeneID),14]
E1 <- EnhancedVolcano(deRes,
                lab = deRes$GeneName,
                x='log2FoldChange',
                y='qvalue',
                title = 'shNorm vs controlNorm',
                subtitle = 'Differential Expression',
                caption = bquote(~Log[2]~ 'foldchange cutoff 0.5; p-value cutoff,0.05'),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                col = c('black','black','green','red'),
                hline = 10e-03,
                labSize = 3,
                pointSize = 3,
                shape = c(1,1,1,25),
                colAlpha = 1,
                legendPosition = 'top',
                legendIconSize = 3,
                drawConnectors = F,
                widthConnectors = 0.5)+
  coord_cartesian(xlim = c(-24,17),ylim = c(0,12.5))
tiff("Plots/shNorm_vs_ct_volcano.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(E1)
dev.off()






#Checking other IHW methods
ihRes <- ihw(pvalue ~ baseMean, data = deRes, alpha = 0.1)
rejections(ihRes)
head(adj_pvalues(ihRes))
sum(adj_pvalues(ihRes) <= 0.1, na.rm = TRUE) == rejections(ihRes)
padjBH <- p.adjust(deRes$pvalue,method = "BH")
sum(padjBH <= 0.1, na.rm = TRUE)
plot(ihRes)
plot(ihRes,what="decisionboundary")

gg <- ggplot(as.data.frame(ihRes), aes(x = pvalue, y = adj_pvalue, col = group)) + 
  geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
gg

gg %+% subset(as.data.frame(ihRes),adj_pvalue <= 0.2)

ihwBonferroni <- ihw(pvalue ~ baseMean, data = deRes, alpha = 0.1, adjustment_type = "bonferroni")


deRes <- na.omit(deRes)
deRes$geneid <- as.numeric(gsub("ENSG[+]*", "", rownames(deRes)))

# set up data frame for ggplotting
rbind(data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$baseMean)/nrow(deRes), 
                 covariate_type="base mean"),
      data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$geneid)/nrow(deRes), 
                 covariate_type="gene id")) %>%
  ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) + 
  facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))


# Plot to show effect of independent filtering
qs <- c(0, quantile(res_shNorm_vs_ctNorm$baseMean[res_shNorm_vs_ctNorm$baseMean > 0], 0:6/6))
bins <- cut(res_shNorm_vs_ctNorm$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(res_shNorm_vs_ctNorm$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")


ggplot(as.data.frame(res_shNorm_vs_ctNorm), aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)