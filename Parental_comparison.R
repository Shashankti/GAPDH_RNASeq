# Comparing parental samples with different glucose levels
# Comparing Norm and Low for control

#PCA plot for norm vs low
pcaData <- plotPCA(vsd[,c(1,2,5,6)],intgroup="condition",returnData=T)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
vst_pca_plot_ln <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")+
  theme_cowplot()
vst_pca_plot_ln



res_pLow_vs_pNorm_2 <-results (dds_2f, contrast = low_parental - norm_parental,independentFiltering = T)
pLow_Norm_df <- as.data.frame(res_pLow_vs_pNorm_2)

#Plot for independent filtering
plot(metadata(res_pLow_vs_pNorm_2)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res_pLow_vs_pNorm_2)$lo.fit, col="red")
abline(v=metadata(res_pLow_vs_pNorm_2)$filterTheta)

# Check the filtering criteria

plot(res_pLow_vs_pNorm_2$baseMean+1, -log10(res_pLow_vs_pNorm_2$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
use <- res_pLow_vs_pNorm_2$baseMean > metadata(res_pLow_vs_pNorm_2)$filterThreshold
h1 <- hist(res_pLow_vs_pNorm_2$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res_pLow_vs_pNorm_2$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

pLow_Norm_df.fixed <- fdrtool::fdrtool(pLow_Norm_df$stat,statistic = "normal")
pLow_Norm_df$qvalue <- p.adjust(pLow_Norm_df.fixed$pval,method = "BH")

#Save csv file as output
write.csv(pLow_Norm_df,file = "Analysis/res_pLow_Norm.csv")

padj.cutoff <- 0.01
lfc.cutoff <- 0.5

threshold <- pLow_Norm_df$qvalue < padj.cutoff & abs(pLow_Norm_df$log2FoldChange) > lfc.cutoff
length(which(threshold))
pLow_Norm_df$threshold <- threshold
pLow_Norm_sig <- data.frame(subset(pLow_Norm_df,threshold==TRUE))
pLow_Norm_sig <- pLow_Norm_sig[order(pLow_Norm_sig$qvalue),]
only_pLow_Norm <- normalized_counts[rownames(pLow_Norm_sig),c(1,2,3,4,5,14)]
only_pLow_Norm$qvalue <- pLow_Norm_sig[match(rownames(only_pLow_Norm),rownames(pLow_Norm_sig)),7]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
only_pLow_Norm_mat <- as.matrix(only_pLow_Norm[,c(2,3,4,5)])
only_pLow_Norm_mat_scaled <- t(scale(t(only_pLow_Norm[,c(2:5)])))
rownames(only_pLow_Norm_mat_scaled) <- only_pLow_Norm$GeneName
ht_pLow <-ComplexHeatmap::Heatmap(na.omit(only_pLow_Norm_mat_scaled),
                                 column_names_side = "bottom",
                                 cluster_rows = T,
                                 cluster_columns = F,
                                 column_names_rot = 45,
                                 column_names_gp = gpar(fontsize=7),
                                 row_names_gp = gpar(fontsize=7),
                                 row_names_rot = 0,
                                 name = "Z-score",width = 4*unit(25,'mm'),
                                 column_title = "controlNorm_vs_controlLow")
ht_pLow
tiff("Plots/ctLow_vs_ctNorm_heatmap",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_pLow)
dev.off()
pLow_Norm_sig$GeneName <- only_pLow_Norm[match(rownames(pLow_Norm_sig),only_pLow_Norm$GeneID),6]

#Save csv file as output
write.csv(pLow_Norm_sig,file = "Analysis/res_pLow_Norm_sig.csv")


# Comparision between parental norm and Gal
res_pGal_vs_pNorm_Unshrunk <-results (dds_2f, contrast = gal_parental - norm_parental,independentFiltering = T)
pGal_Norm_df <- as.data.frame(res_pGal_vs_pNorm_Unshrunk)

#Plot for independent filtering
plot(metadata(res_pGal_vs_pNorm_Unshrunk)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res_pGal_vs_pNorm_Unshrunk)$lo.fit, col="red")
abline(v=metadata(res_pGal_vs_pNorm_Unshrunk)$filterTheta)

# Check the filtering criteria

plot(res_pGal_vs_pNorm_Unshrunk$baseMean+1, -log10(res_pGal_vs_pNorm_Unshrunk$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
use <- res_pGal_vs_pNorm_Unshrunk$baseMean > metadata(res_pGal_vs_pNorm_Unshrunk)$filterThreshold
h1 <- hist(res_pGal_vs_pNorm_Unshrunk$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res_pGal_vs_pNorm_Unshrunk$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

pGal_Norm_df.fixed <- fdrtool::fdrtool(pGal_Norm_df$stat,statistic = "normal")
pGal_Norm_df$qvalue <- p.adjust(pGal_Norm_df.fixed$pval,method = "BH")

#Save csv file as output
write.csv(pGal_Norm_df,file = "Analysis/res_pGal_Norm.csv")

padj.cutoff <- 0.01
lfc.cutoff <- 0.5

threshold <- pGal_Norm_df$qvalue < padj.cutoff & abs(pGal_Norm_df$log2FoldChange) > lfc.cutoff
length(which(threshold))
pGal_Norm_df$threshold <- threshold
pGal_Norm_sig <- data.frame(subset(pGal_Norm_df,threshold==TRUE))
pGal_Norm_sig <- pGal_Norm_sig[order(pGal_Norm_sig$qvalue),]
only_pGal_Norm <- normalized_counts[rownames(pGal_Norm_sig),c(1,2,3,6,7,14)]
only_pGal_Norm$qvalue <- pGal_Norm_sig[match(rownames(only_pGal_Norm),rownames(pGal_Norm_sig)),7]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
only_pGal_Norm_mat <- as.matrix(only_pGal_Norm[,c(2,3,4,5)])
only_pGal_Norm_mat_scaled <- t(scale(t(only_pGal_Norm[,c(2:5)])))
rownames(only_pGal_Norm_mat_scaled) <- only_pGal_Norm$GeneName
ht_pGal <-ComplexHeatmap::Heatmap(na.omit(only_pGal_Norm_mat_scaled),
                                 column_names_side = "bottom",
                                 cluster_rows = T,
                                 cluster_columns = F,
                                 column_names_rot = 45,
                                 column_names_gp = gpar(fontsize=7),
                                 row_names_gp = gpar(fontsize=7),
                                 row_names_rot = 0,
                                 name = "Z-score",width = 4*unit(25,'mm'),
                                 column_title = "controlNorm_vs_controlGal")
ht_pGal
tiff("Plots/ctGal_vs_ctNorm_heatmap",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_pGal)
dev.off()

pGal_Norm_sig$GeneName <- only_pGal_Norm[match(rownames(pGal_Norm_sig),only_pGal_Norm$GeneID),6]

#Save csv file as output
write.csv(pGal_Norm_sig,file = "Analysis/res_pGal_Norm_sig.csv")