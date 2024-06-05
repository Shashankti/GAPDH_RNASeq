# Running multiple testing correction for Gal and Low

library(dplyr)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(EnhancedVolcano)

#Import deSeq object at df
res_shLow_vs_ctLow_unshrunk <- results(dds_2f, contrast = low_shRNA - low_parental,format = "DataFrame")
res_shGal_vs_ctGal_unshrunk <- results(dds_2f, contrast = gal_shRNA - gal_parental)
#convert to dataframe
Low_df <- as.data.frame(res_shLow_vs_ctLow_unshrunk)
Gal_df <- as.data.frame(res_shGal_vs_ctGal_unshrunk)

#Apply fdrtool to value
Low_df.fixed <- fdrtool::fdrtool(Low_df$stat,statistic = "normal")
Low_df$qvalue <- p.adjust(Low_df.fixed$pval,method = "BH")

Gal_df.fixed <- fdrtool::fdrtool(Gal_df$stat,statistic = "normal")
Gal_df$qvalue <- p.adjust(Gal_df.fixed$pval,method = "BH")

#save csv output
write.csv(Low_df,file = "res_shLow_vs_ct_Low")
write.csv(Gal_df,file = "res_shGal_vs_ct_Gal")

#Making heatmaps for Low Glucose
padj.cutoff <- 0.05
lfc.cutoff <- 0.5

threshold <- Low_df$qvalue < padj.cutoff & abs(Low_df$log2FoldChange) > lfc.cutoff
length(which(threshold))
Low_df$threshold <- threshold
Low_sig <- data.frame(subset(Low_df,threshold==TRUE))
Low_sig <- Low_sig[order(Low_sig$qvalue),]
only_low <- normalized_counts[rownames(Low_sig),c(1,4,5,9,12,14)]
only_low$qvalue <- Low_sig[match(rownames(only_low),rownames(Low_sig)),7]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
only_low_mat <- as.matrix(only_low[,c(2,3,4,5)])
only_low_mat_scaled <- t(scale(t(only_low[,c(2:5)])))
rownames(only_low_mat_scaled) <- only_low$GeneName
Low_sig$GeneName <- only_low[match(rownames(Low_sig),only_low$GeneID),6]

#save output as csv file
write.csv(Low_sig,file = "res_shLow_vs_ct_Low_sig")


ht_Low <-ComplexHeatmap::Heatmap(na.omit(only_low_mat_scaled),
                                  column_names_side = "bottom",
                                  cluster_rows = T,
                                  cluster_columns = F,
                                  column_names_rot = 45,
                                  column_names_gp = gpar(fontsize=7),
                                  row_names_gp = gpar(fontsize=5),
                                 row_names_rot = 0,
                                 name = "Z-score",width = 4*unit(25,'mm'),
                                 column_title = "controlLow_vs_shLow")
ht_Low
tiff("Plots/shRNA_vs_ctLow_heatmap",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_Low)
dev.off()
# Making plots for Gal comparision
padj.cutoff <- 0.05
lfc.cutoff <- 0.5


#making a volcano plot
Low_df$GeneName <- normalized_counts[match(rownames(Low_df),normalized_counts$GeneID),14]
E1 <- EnhancedVolcano(Low_df,
                      lab = Low_df$GeneName,
                      x='log2FoldChange',
                      y='qvalue',
                      title = 'shLow vs controlLow',
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
tiff("Plots/shLow_vs_ct_volcano.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(E1)
dev.off()



threshold <- Gal_df$qvalue < padj.cutoff & abs(Gal_df$log2FoldChange) > lfc.cutoff
length(which(threshold))
Gal_df$threshold <- threshold
Gal_sig <- data.frame(subset(Gal_df,threshold==TRUE))
Gal_sig <- Gal_sig[order(Gal_sig$qvalue),]
only_Gal <- normalized_counts[rownames(Gal_sig),c(1,6,7,10,13,14)]
only_Gal$qvalue <- Gal_sig[match(rownames(only_Gal),rownames(Gal_sig)),7]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
only_Gal_mat <- as.matrix(only_Gal[,c(2,3,4,5)])
only_Gal_mat_scaled <- t(scale(t(only_Gal[,c(2:5)])))
rownames(only_Gal_mat_scaled) <- only_Gal$GeneName
Gal_sig$GeneName <- only_Gal[match(rownames(Gal_sig),only_Gal$GeneID),6]

#save output as csv file
write.csv(Gal_sig,file = "res_shLow_vs_ct_Gal_sig")

ht_Gal <-ComplexHeatmap::Heatmap(na.omit(only_Gal_mat_scaled),
                                 column_names_side = "bottom",
                                 cluster_rows = T,
                                 cluster_columns = F,
                                 column_names_rot = 45,
                                 column_names_gp = gpar(fontsize=7),
                                 row_names_gp = gpar(fontsize=5),
                                 row_names_rot = 0,
                                 name = "Z-score",width = 4*unit(25,'mm'),
                                 column_title = "controlGal_vs_shGal")
ht_Gal
tiff("Plots/shRNA_vs_ctGal_heatmap",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_Gal)
dev.off()

#making volcano plot
Gal_df$GeneName <- normalized_counts[match(rownames(Gal_df),normalized_counts$GeneID),14]
E1 <- EnhancedVolcano(Gal_df,
                      lab = Gal_df$GeneName,
                      x='log2FoldChange',
                      y='qvalue',
                      title = 'shGal vs controlGal',
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
  coord_cartesian(xlim = c(-20,20),ylim = c(0,12.5))
tiff("Plots/shGal_vs_ct_volcano.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(E1)
dev.off()


# Comparison between sh all conditions and control
res_control_vs_sh <- results(dds_2f,c("type","shRNA","parental"))
sh_vs_control_df <- as.data.frame(res_control_vs_sh)


#save the csv file
write.csv(sh_vs_control_df,file = "res_shRNA_vs_control")


# distribution of p-value
ggplot(as.data.frame(res_control_vs_sh), aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0,fill="azure4")+
  theme_cowplot()+
  ggtitle("Control vs shRNA p-value distribution")
#Info regarding metadata
metadata(res_control_vs_sh)$filterThreshold

#Plot for independent filtering
plot(metadata(res_control_vs_sh)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res_control_vs_sh)$lo.fit, col="red")
abline(v=metadata(res_control_vs_sh)$filterTheta)

# Check the filtering criteria

plot(res_control_vs_sh$baseMean+1, -log10(res_control_vs_sh$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))
use <- res_control_vs_sh$baseMean > metadata(res_control_vs_sh)$filterThreshold
h1 <- hist(res_control_vs_sh$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res_control_vs_sh$pvalue[use], breaks=0:50/50, plot=FALSE)

colori <- c(`do not pass`="khaki", `pass`="powderblue")

barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

# since the counts are sufficiently high we do not need to apply the fdr approach to this data
# however we still try to use the fdr apporach for consistency 

padj.cutoff <- 0.1
lfc.cutoff <- 0.5

threshold <- sh_vs_control_df$padj < padj.cutoff & abs(sh_vs_control_df$log2FoldChange) > lfc.cutoff
length(which(threshold))

sh_vs_control_df$threshold <- threshold
sh_vs_control_sig <- data.frame(subset(sh_vs_control_df,threshold==TRUE))
sh_vs_control_sig <- sh_vs_control_sig[order(sh_vs_control_sig$padj),]
all_counts <- normalized_counts[rownames(sh_vs_control_sig),]
all_counts$padj <- sh_vs_control_sig[match(rownames(all_counts),rownames(sh_vs_control_sig)),6]
#melted_norm <- data.frame(reshape::melt(sigOE_ordered_qValue)) 
all_counts_mat <- as.matrix(only_Gal[,c(2,3,4,5)])
all_counts_mat_scaled <- t(scale(t(all_counts[,c(2:13)])))
rownames(all_counts_mat_scaled) <- all_counts$GeneName

sh_vs_control_sig$GeneName <- all_counts[match(rownames(sh_vs_control_sig),all_counts$GeneID),6]

#save output as csv file
write.csv(sh_vs_control_sig,file = "res_shRNA_vs_ct_sig")

ht_ct_sh <-ComplexHeatmap::Heatmap(all_counts_mat_scaled,
                                 column_names_side = "bottom",
                                 cluster_rows = T,
                                 cluster_columns = F,
                                 column_names_rot = 45,
                                 column_names_gp = gpar(fontsize=7),
                                 row_names_gp = gpar(fontsize=5),
                                 row_names_rot = 0,
                                 name = "Z-score",
                                 column_title = "control_vs_shRNA",
                                 column_split = meta$sampletype)
ht_ct_sh
tiff("Plots/shRNA_vs_ct_heatmap",width = 35,height = 21,units = 'cm',res = 300)
plot(ht_ct_sh)
dev.off()

# making volcano plots

sh_vs_control_df$GeneName <- normalized_counts[match(rownames(sh_vs_control_df),normalized_counts$GeneID),14]
E1 <- EnhancedVolcano(sh_vs_control_df,
                      lab = sh_vs_control_df$GeneName,
                      x='log2FoldChange',
                      y='padj',
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
tiff("Plots/shLow_vs_ct_volcano.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(E1)
dev.off()
