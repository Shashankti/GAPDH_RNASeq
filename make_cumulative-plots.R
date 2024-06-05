#!mkae cumulative distribution plots


library(ggplot2)
library(cowplot)
library(ggrepel)
library(data.table)
library(rstatix)
library(plyr)
library(tidyverse)
library(dbplyr)
library(ggpubr)

`%nin%` = Negate(`%in%`)
#load in the clip peaks file



# old peaks,
#' Ignore in future
Gal_joined <- fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Union/Gal_joined_annotated.bed")
Norm_joined <- fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Union/Norm_joined_annotated.bed")
Lo_joined <- fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Union/Lo_joined_annotated.bed")


#new peaks joined
#Use these peaks for future

Gal_joined <- fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Union/JoinedPeaks/Joined_Gal_annotated.bed")
Norm_joined <- fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Union/JoinedPeaks/Joined_Norm_annotated.bed")
Lo_joined <- fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Union/JoinedPeaks/Joined_Lo_annotated.bed")

#-----------------------------------------------------------#

Gal_ovrlap <- fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Merged_Gal_annotated.bed")
Norm_ovrlap <- fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Merged_Norm_annotated.bed")
Lo_ovrlap <- fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Merged_Lo_annotated.bed")
#make genes_list
clip_gal_genes <- unique(Gal_joined$V4)
clip_lo_genes <- unique(Lo_joined$V4)
clip_norm_genes <- unique(Norm_joined$V4)
#make genes_list
clip_gal_genes <- unique(Gal_ovrlap$V5)
clip_lo_genes <- unique(Lo_ovrlap$V5)
clip_norm_genes <- unique(Norm_ovrlap$V5)

#make cdf for all genees
#reads the rlog norm count for CDF
# rld <- assay(rlog(dds_2f))
#make a combined dataframe

#Apply fdrtool to value
#we rerun the multiple testing correction to remove sparse values with NAs in 
df_res_control_vs_shRNA <- as.data.frame(res_control_vs_sh)
df_res_control_vs_shRNA.fixed <- fdrtool::fdrtool(df_res_control_vs_shRNA$stat,statistic = "normal")
df_res_control_vs_shRNA$qvalue <- p.adjust(df_res_control_vs_shRNA.fixed$pval,method = "BH")
df_res_control_vs_shRNA$GeneName <- gene_IDs[match(rownames(df_res_control_vs_shRNA),gene_IDs[,1]),2]


#make new volcano plot
E1 <- EnhancedVolcano(df_res_control_vs_shRNA,
                      lab = df_res_control_vs_shRNA$GeneName,
                      x='log2FoldChange',
                      y='qvalue',
                      title = 'shRNA vs control',
                      subtitle ='',
                      subtitleLabSize = 0.1,
                      caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.1'),
                      pCutoff = 0.1,
                      FCcutoff = 1,
                      labSize = 3,
                      pointSize = 0.7,
                      colAlpha = 0.8,
                      legendPosition = 'top',
                      legendIconSize = 2,
                      drawConnectors = F,
                      widthConnectors = 0.5)+
  coord_cartesian(xlim = c(-16,16),
                  ylim = c(0,15),
                  clip = 'off')
tiff("Plots/Control_vs_shRNA_volcano.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(E1)
dev.off()
write.table(df_res_control_vs_shRNA,"all_control_vs_shRNA_ko.tsv",sep = "\t",quote = F)

deRes <- arrange(d,qvalue) %>% rownames_to_column("GeneID")

p.adj <- p.adjust(res_shNorm_vs_ctNorm_Unshrunk$pvalue,method = "BH")

## adding negative control with random genes 
n <- dim(df_res_control_vs_shRNA)[1]
sampl <- sample.int(n, 400)

df_sampld <- df_res_control_vs_shRNA[sampl,]
df_sampld[rownames(df_sampld) %nin% unique(union(union(clip_gal_genes,clip_lo_genes),clip_norm_genes)),]

logFoldChanges <- c(df_res_control_vs_shRNA[df_res_control_vs_shRNA$qvalue<10,]$log2FoldChange,
                    res_control_vs_sh[rownames(res_control_vs_sh)%in% clip_gal_genes,2],
                    res_control_vs_sh[rownames(res_control_vs_sh)%in% clip_norm_genes,2],
                    res_control_vs_sh[rownames(res_control_vs_sh)%in% clip_lo_genes,2],
                    df_sampld$log2FoldChange)

types <- c(rep("All",dim(df_res_control_vs_shRNA[df_res_control_vs_shRNA$qvalue<10,])[1],),
           rep("Gal", length(res_control_vs_sh[rownames(res_control_vs_sh)%in% clip_gal_genes,2])),
           rep("Norm",length(res_control_vs_sh[rownames(res_control_vs_sh)%in% clip_norm_genes,2])),
           rep("Low",length(res_control_vs_sh[rownames(res_control_vs_sh)%in% clip_lo_genes,2])),
           rep("NonTargets",length(df_sampld$log2FoldChange)))
           
comb_df <- data.frame(types,logFoldChanges)
comb_df_summary <- ddply(comb_df, .(types), summarise, median = median(logFoldChanges))

pvalue_list = NULL
pvalue_list[1] <- 1
pvalue_list[4] <- wilcox.test(comb_df[comb_df$types=="All",2], comb_df[comb_df$types=="NonTargets",2])$p.value
comb_df_summary <- cbind(comb_df_summary, pvalue_list)
comb_df_summary

fold_change_list = NULL
fold_change_list[1] <- NA

fold_change_list[2] <- mean(comb_df[comb_df$types=="Gal",]$logFoldChanges) - mean(comb_df[comb_df$types=="All",]$logFoldChanges)
fold_change_list[3] <- mean(comb_df[comb_df$types=="Low",]$logFoldChanges) - mean(comb_df[comb_df$types=="All",]$logFoldChanges)
fold_change_list[4] <- mean(comb_df[comb_df$types=="NonTargets",]$logFoldChanges) - mean(comb_df[comb_df$types=="All",]$logFoldChanges)
fold_change_list[5] <- mean(comb_df[comb_df$types=="Norm",]$logFoldChanges) - mean(comb_df[comb_df$types=="All",]$logFoldChanges)

comb_df_summary <- cbind(comb_df_summary,fold_change_list)
comb_df_summary
# attach(loadNamespace("enrichplot"), name = "enrichplot_all")


#set up a color pallette
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)

cols <- c("All"="#3B3B3BFF","Gal"="#d35400","Low"="#2980b9","Norm"="#27ae60","NonTargets"="#8F7700FF")

tiff("Plots/Cumulative_distribution/JoinedPeaks_cumulative_dist_with_all.tiff",width = 35,height = 21,units = 'cm',res = 300)
ggplot(comb_df,aes(logFoldChanges,colour = types)) +
  stat_ecdf(size=0.8)+
  scale_color_manual(values = cols) + 
  theme_cowplot()+
  ggtitle("CDF for joined peaks with all genes")+
  xlim(-2,2) + ylab("Cumulative Distribution")+
  labs(color="Type")+
  annotate("text", x= -1.8, y= 0.97, label = "21,504 all genes,",color = "#3B3B3BFF")+
  annotate("text", x= -1.8, y= 0.95, label = "Median=0.006",color = "#3B3B3BFF")+
  annotate("text",x=-1.8,y=0.89, label = "324 Gal Genes,",color = "#d35400") + 
  annotate("text",x=-1.8,y=0.87, label = "Median=-0.058",color = "#d35400") + 
  annotate("text",x=-1.8,y=0.85, label = "(FC=-0.2394,", color = "#d35400") + 
  annotate("text",x=-1.8,y=0.83, label = "P = 6.44e-04)", color = "#d35400") + 
  annotate("text",x=-1.8,y=0.77, label = "449 Low Genes", color = "#2980b9") +
  annotate("text",x=-1.8,y=0.75, label = "Median= -0.093", color = "#2980b9") +
  annotate("text",x=-1.8,y=0.73,label = "(FC=-0.2282,", color = "#2980b9") +
  annotate("text",x=-1.8,y=0.71,label = "P=9.35e-06)",color = "#2980b9") +
  annotate("text",x=-1.8,y=0.67,label= "375 Norm Genes", color = "#27ae60") +
  annotate("text",x=-1.8,y=0.65,label= "Median=-0.2643", color = "#27ae60") +
  annotate("text",x=-1.8,y=0.63,label = "(FC=-0.2585118,", color = "#27ae60") + 
  annotate("text",x=-1.8,y=0.61,label = "P=5.043-06)", color = "#27ae60")+
  annotate("text",x=-1.8,y=0.57,label= "Non Targets", color = "#8F7700FF") +
  annotate("text",x=-1.8,y=0.55,label= "Median=0.011", color = "#8F7700FF") +
  annotate("text",x=-1.8,y=0.53,label = "(FC=0.0404,", color = "#8F7700FF") + 
  annotate("text",x=-1.8,y=0.51,label = "P=0.66)", color = "#8F7700FF")
dev.off()


#only Gal plot


cols <- c("All"="#3B3B3BFF","Gal"="#d35400","NegControl"="#8F7700FF")
tiff("Plots/Cumulative_distribution/JoinedPeaks_cumulative_dist_only_gal.tiff",width = 35,height = 21,units = 'cm',res = 300)
ggplot(comb_df[comb_df$types %in% c("All","Gal","NegControl"),],aes(logFoldChanges,colour = types)) +
  stat_ecdf(size=0.8)+
  scale_color_manual(values = cols) + 
  theme_cowplot()+
  ggtitle("CDF for joined peaks with Gal Peaks")+
  xlim(-2,2) + ylab("Cumulative Distribution")+
  labs(color="Type")+
  annotate("text", x= -1.8, y= 0.97, label = "21,504 all genes,",color = "#3B3B3BFF")+
  annotate("text", x= -1.8, y= 0.95, label = "Median=0.006",color = "#3B3B3BFF")+
  annotate("text",x=-1.8,y=0.89, label = "324 Gal Genes,",color = "#d35400") + 
  annotate("text",x=-1.8,y=0.87, label = "Median=-0.058",color = "#d35400") + 
  annotate("text",x=-1.8,y=0.85, label = "(FC=-0.2394,", color = "#d35400") + 
  annotate("text",x=-1.8,y=0.83, label = "P = 6.44e-04)", color = "#d35400") + 
  annotate("text",x=-1.8,y=0.77,label= "Negative Control", color = "#8F7700FF") +
  annotate("text",x=-1.8,y=0.75,label= "Median=0.011", color = "#8F7700FF") +
  annotate("text",x=-1.8,y=0.73,label = "(FC=0.0404,", color = "#8F7700FF") + 
  annotate("text",x=-1.8,y=0.71,label = "P=0.66)", color = "#8F7700FF")
dev.off()


#---------------------------------------------------------------------------------------#
#only low clip targets 
cols <- c("All"="#3B3B3BFF","Lo"="#2980b9","NegControl"="#8F7700FF")

tiff("Plots/Cumulative_distribution/JoinedPeaks_cumulative_dist_only_low.tiff",width = 35,height = 21,units = 'cm',res = 300)
ggplot(comb_df[comb_df$types %in% c("All","Lo","NegControl"),],aes(logFoldChanges,colour = types)) +
  stat_ecdf(size=0.8)+
  scale_color_manual(values = cols) + 
  theme_cowplot()+
  ggtitle("CDF for joined peaks with Low Glucose Peaks")+
  xlim(-2,2) + ylab("Cumulative Distribution")+
  labs(color="Type")+
  annotate("text", x= -1.8, y= 0.97, label = "21,504 all genes,",color = "#3B3B3BFF")+
  annotate("text", x= -1.8, y= 0.95, label = "Median=0.006",color = "#3B3B3BFF")+
  annotate("text",x=-1.8,y=0.87, label = "449 Low Genes", color = "#2980b9") +
  annotate("text",x=-1.8,y=0.85, label = "Median= -0.093", color = "#2980b9") +
  annotate("text",x=-1.8,y=0.83,label = "(FC=-0.2282,", color = "#2980b9") +
  annotate("text",x=-1.8,y=0.81,label = "P=9.35e-06)",color = "#2980b9") +
  annotate("text",x=-1.8,y=0.77,label= "Negative Control", color = "#8F7700FF") +
  annotate("text",x=-1.8,y=0.75,label= "Median=0.011", color = "#8F7700FF") +
  annotate("text",x=-1.8,y=0.73,label = "(FC=0.0404,", color = "#8F7700FF") + 
  annotate("text",x=-1.8,y=0.71,label = "P=0.66)", color = "#8F7700FF")
dev.off()


#------------------------------------------------------------------------------------------------------------------#
#only norm glucose clip peaks

cols <- c("All"="#3B3B3BFF","Norm"="#27ae60","NegControl"="#8F7700FF")

tiff("Plots/Cumulative_distribution/JoinedPeaks_cumulative_dist_only_norm.tiff",width = 35,height = 21,units = 'cm',res = 300)
ggplot(comb_df[comb_df$types %in% c("All","Norm","NegControl"),],aes(logFoldChanges,colour = types)) +
  stat_ecdf(size=0.8)+
  scale_color_manual(values = cols) + 
  theme_cowplot()+
  ggtitle("CDF for joined peaks with Normal Glucose Peaks")+
  xlim(-2,2) + ylab("Cumulative Distribution")+
  labs(color="Type")+
  annotate("text", x= -1.8, y= 0.97, label = "21,504 all genes,",color = "#3B3B3BFF")+
  annotate("text", x= -1.8, y= 0.95, label = "Median=0.006",color = "#3B3B3BFF")+
  annotate("text",x=-1.8,y=0.87,label= "375 Norm Genes", color = "#27ae60") +
  annotate("text",x=-1.8,y=0.85,label= "Median=-0.2643", color = "#27ae60") +
  annotate("text",x=-1.8,y=0.83,label = "(FC=-0.2585118,", color = "#27ae60") + 
  annotate("text",x=-1.8,y=0.81,label = "P=5.043-06)", color = "#27ae60")+
  annotate("text",x=-1.8,y=0.77,label= "Negative Control", color = "#8F7700FF") +
  annotate("text",x=-1.8,y=0.75,label= "Median=0.011", color = "#8F7700FF") +
  annotate("text",x=-1.8,y=0.73,label = "(FC=0.0404,", color = "#8F7700FF") + 
  annotate("text",x=-1.8,y=0.71,label = "P=0.66)", color = "#8F7700FF")
dev.off()


#-------------------------------------------------------------------------------------------------------------------#
#only Gal vs Gal

gal_fc <- c(Gal_df$log2FoldChange, Gal_df[rownames(Gal_df) %in% clip_gal_genes,2])
types_gal <- c(rep("All",dim(Gal_df)[1]), rep("Clip",length(Gal_df[rownames(Gal_df) %in% clip_gal_genes,2])))
comb_gal <- data.frame(types_gal, gal_fc)
comb_gal_summary <- ddply(comb_gal, .(types_gal), summarise, median = median(gal_fc))
pvalue_gal <- wilcox.test(comb_gal[comb_gal$types_gal=="All",2], comb_gal[comb_gal$types_gal=="Clip",2])$p.value


comb_df %>% 
  group_by(types) %>%
  get_summary_stats(logFoldChanges, type = "mean_sd")
stat.test <- comb_df %>%
  cohens_d(logFoldChanges ~ types) %>%
  add_significance()


png("~/DEseq2_RNAseq_GAPDH_KD/only_gal_cdf_merged.png",width = 1200,height = 800)
ggplot(comb_gal,aes(gal_fc,colour = types_gal)) +
  stat_ecdf()+
  scale_color_hue(name="Types",labels = c('All','Clip')) + 
  theme_cowplot()+
  ggtitle("CDF for common peaks in Gal ")+
  xlim(-2,2) + ylab("Cumulative Distribution")+
  geom_label_repel(data=comb_gal_summary, aes(x = median, label=median, y= 0.05),
                   size = 3, segment.size=0.25, nudge_x = 0.5, direction = "y",
                   hjust=0) + geom_point(data = comb_gal_summary,aes(x = median, y= 0.05, shape = types_gal), size=2)
dev.off()


#only Low
lo_fc <- c(Low_df$log2FoldChange, Low_df[rownames(Low_df) %in% clip_lo_genes,2])
types_lo <- c(rep("All",dim(Low_df)[1]), rep("Clip",length(Low_df[rownames(Low_df) %in% clip_lo_genes,2])))
comb_lo <- data.frame(types_lo, lo_fc)
comb_lo_summary <- ddply(comb_lo, .(types_lo), summarise, median = median(lo_fc))
pvalue_lo <- wilcox.test(comb_lo[comb_lo$types_lo=="All",2], comb_lo[comb_lo$types_lo=="Clip",2])$p.value



png("~/DEseq2_RNAseq_GAPDH_KD/only_lo_cdf_merged.png",width = 1200,height = 800)
ggplot(comb_lo,aes(lo_fc,colour = types_lo)) +
  stat_ecdf()+
  scale_color_hue(name="Types",labels = c('All','Clip')) + 
  theme_cowplot()+
  ggtitle("CDF for common peaks in Low ")+
  xlim(-2,2) + ylab("Cumulative Distribution")+
  geom_label_repel(data=comb_lo_summary, aes(x = median, label=median, y= 0.05),
                   size = 3, segment.size=0.25, nudge_x = 0.5, direction = "y",
                   hjust=0) + geom_point(data = comb_lo_summary,aes(x = median, y= 0.05, shape = types_lo), size=2)
dev.off()



#3#only norm
Norm_df <- deRes
norm_fc <- c(Norm_df$log2FoldChange, Norm_df[rownames(Low_df) %in% clip_norm_genes,2])
types_norm <- c(rep("All",dim(Norm_df)[1]), rep("Clip",length(Norm_df[rownames(Norm_df) %in% clip_norm_genes,2])))
comb_norm <- data.frame(types_norm, norm_fc)
comb_norm_summary <- ddply(comb_norm, .(types_norm), summarise, median = median(norm_fc))
pvalue_norm <- wilcox.test(comb_norm[comb_norm$types_norm=="All",2], comb_norm[comb_norm$types_norm=="Clip",2])$p.value



png("~/DEseq2_RNAseq_GAPDH_KD/only_norm_cdf_merged.png",width = 1200,height = 800)
ggplot(comb_norm,aes(norm_fc,colour = types_norm)) +
  stat_ecdf()+
  scale_color_hue(name="Types",labels = c('All','Clip')) + 
  theme_cowplot()+
  ggtitle("CDF for common peaks in Norm ")+
  xlim(-2,2) + ylab("Cumulative Distribution")+
  geom_label_repel(data=comb_norm_summary, aes(x = median, label=median, y= 0.05),
                   size = 3, segment.size=0.25, nudge_x = 0.5, direction = "y",
                   hjust=0) + geom_point(data = comb_norm_summary,aes(x = median, y= 0.05, shape = types_norm), size=2)
dev.off()





#--------------------------------------------------------------------------------------------------------------------
#make cumulative distribution plot with only the consensus peaks for clarity
all_joined <- fread("~/CLIP_Seq/Data_Strict_Align/Peaks/PureClip/Union/JoinedPeaks/All_int_annotated.bed")
#make genes_list
clip_genes <- unique(all_joined$V4)


all_logFoldChanges <- c(df_res_control_vs_shRNA[df_res_control_vs_shRNA$qvalue<10,]$log2FoldChange,
                        res_control_vs_sh[rownames(res_control_vs_sh)%in% clip_genes,2],
                        df_sampld$log2FoldChange)

all_types <- c(rep("All",dim(df_res_control_vs_shRNA[df_res_control_vs_shRNA$qvalue<10,])[1],),
               rep("Consensus", length(res_control_vs_sh[rownames(res_control_vs_sh)%in% clip_genes,2])),
               rep("NegControl",length(df_sampld$log2FoldChange)))

all_df <- data.frame(all_types,all_logFoldChanges)
all_df_summary <- ddply(all_df, .(all_types), summarise, median = median(all_logFoldChanges))

pvalue_list = NULL
pvalue_list[1] <- 1
pvalue_list[3] <- wilcox.test(all_df[all_df$all_types=="All",2], all_df[all_df$all_types=="NegControl",2])$p.value
pvalue_list[2] <- wilcox.test(all_df[all_df$all_types=="All",2], all_df[all_df$all_types=="Consensus",2])$p.value

all_df_summary <- cbind(all_df_summary, pvalue_list)
all_df_summary

fold_change_list = NULL
fold_change_list[1] <- NA

fold_change_list[2] <- mean(all_df[all_df$all_types=="Consensus",]$all_logFoldChanges) - 
  mean(all_df[all_df$all_types=="All",]$all_logFoldChanges)
fold_change_list[3] <- mean(all_df[all_df$all_types=="NegControl",]$all_logFoldChanges) - 
  mean(all_df[all_df$all_types=="All",]$all_logFoldChanges)


all_df_summary <- cbind(all_df_summary,fold_change_list)
all_df_summary



tiff("Plots/Cumulative_distribution/JoinedPeaks_cumulative_dist_with_consensus.tiff",width = 35,height = 21,units = 'cm',res = 300)
ggplot(all_df,aes(all_logFoldChanges,colour = all_types)) +
  stat_ecdf(size=0.8)+
  scale_color_manual(values = cols) + 
  theme_cowplot()+
  ggtitle("CDF for joined peaks with all genes")+
  xlim(-2,2) + ylab("Cumulative Distribution")+
  labs(color="Type")+
  annotate("text", x= -1.8, y= 0.97, label = "21,504 all genes,",color = "#3B3B3BFF")+
  annotate("text", x= -1.8, y= 0.95, label = "Median=0.006",color = "#3B3B3BFF")+
  annotate("text",x=-1.8,y=0.89, label = "324 Gal Genes,",color = "#d35400") + 
  annotate("text",x=-1.8,y=0.87, label = "Median=-0.058",color = "#d35400") + 
  annotate("text",x=-1.8,y=0.85, label = "(FC=-0.2394,", color = "#d35400") + 
  annotate("text",x=-1.8,y=0.83, label = "P = 6.44e-04)", color = "#d35400") + 
  annotate("text",x=-1.8,y=0.77, label = "449 Low Genes", color = "#2980b9") +
  annotate("text",x=-1.8,y=0.75, label = "Median= -0.093", color = "#2980b9") +
  annotate("text",x=-1.8,y=0.73,label = "(FC=-0.2282,", color = "#2980b9") +
  annotate("text",x=-1.8,y=0.71,label = "P=9.35e-06)",color = "#2980b9") +
  annotate("text",x=-1.8,y=0.67,label= "375 Norm Genes", color = "#27ae60") +
  annotate("text",x=-1.8,y=0.65,label= "Median=-0.2643", color = "#27ae60") +
  annotate("text",x=-1.8,y=0.63,label = "(FC=-0.2585118,", color = "#27ae60") + 
  annotate("text",x=-1.8,y=0.61,label = "P=5.043-06)", color = "#27ae60")+
  annotate("text",x=-1.8,y=0.57,label= "Negative Control", color = "#8F7700FF") +
  annotate("text",x=-1.8,y=0.55,label= "Median=0.011", color = "#8F7700FF") +
  annotate("text",x=-1.8,y=0.53,label = "(FC=0.0404,", color = "#8F7700FF") + 
  annotate("text",x=-1.8,y=0.51,label = "P=0.66)", color = "#8F7700FF")
dev.off()




# detachAllPackages <- function() {
#   
#   basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
#   
#   package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
#   
#   package.list <- setdiff(package.list,basic.packages)
#   
#   if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
#   
# }
# 
# detachAllPackages()
# 
# 
# unload_pkgs <- \(exc=NULL) {
#   bpk <- c("compiler", "graphics", "tools", "utils", "grDevices", "stats", 
#            "datasets", "methods", "base", "rstudioapi") |> c(exc)
#   while (length(setdiff(loadedNamespaces(), bpk)) > 0) {
#     lapply(setdiff(loadedNamespaces(), bpk), \(x) {
#       try(unloadNamespace(x), silent=TRUE)
#     })
#   }
# }
# unload_pkgs()