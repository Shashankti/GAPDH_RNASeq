#correlation between TCGA GAPDH vs GAPDH-CLIP targets



#import the counts file

library(data.table)
library(rstatix)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(cowplot)
library(DGEobj.utils)
library(ggpubr)

#get gene lenghts
library(GenomicFeatures)
txdb <- GenomicFeatures::makeTxDbFromGFF("~/CLIP_Seq/Data/Genome/Homo_sapiens.GRCh38.108.gtf",format = "gtf")
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
df <- as.data.frame(exonic.gene.sizes[rownames(TCGA_df),])
rownames(df) <- rownames(TCGA_df)
#tpm counts 
#TCGA LAML raw counts
TCGA_df
#generate TPM normalized counts
nm_tcga <- convertCounts(TCGA_df,
                         unit = "TPM",
                         geneLength = as.matrix(df),
                         log = FALSE
                         )
nm_tcga <- as.data.frame(nm_tcga)
nm_tcga$GeneName <-  gene_IDs[match(rownames(nm_tcga), gene_IDs[,1]),2]

imp <- nm_tcga[nm_tcga$GeneName %in% c("GAPDH", "PKM"),]
imp <- as.data.frame(t(imp))
colnames(imp) <- c("PKM", "GAPDH")
imp <- head(imp,-1)
imp[,1] <- as.numeric(imp[,1])
imp[,2] <- as.numeric(imp[,2])
#calculate correlation
cor.test(imp$PKM,imp$GAPDH,method = "spearman")



p1 <- ggplot(imp, aes(x=GAPDH,y=PKM))+
  geom_point(size=2,alpha=0.5)+
  theme_cowplot()+
  ylab("PKM Normalized Counts")+
  ggtitle("GAPDH vs. PKM")+
  theme(plot.title = element_text(size = 15),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.x = element_text(margin = margin(t(30))))+
  xlab(expression(paste("         GAPDH Normalized counts \n Spearmann p<2.2e-16 coeff=0.60(n=151)")))
coord_cartesian(ylim = c(0,30000), clip = "off")
tiff("GAPDH_vs_PKM.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(p1)
dev.off()





#scatter plot
sp <- ggscatter(imp,x="GAPDH",y="RPL13A",
                add = "reg.line",
                add.params = list(color="blue",fill="lightgray"),
                conf.int = TRUE,
                size=2,alpha=0.5)+
  theme_cowplot()+
  ylab("RPL13A Normalized Counts")+
  ggtitle("GAPDH vs. RPL13A")+
  theme(plot.title = element_text(size = 15),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        axis.title.x = element_text(margin = margin(t(30))))

sp + stat_cor(method = "spearman")

ggplot(imp, aes(x = ,y = log10(Normalized), fill = Type))+
  geom_violin(trim = FALSE, width=0.4)+
  theme_cowplot()+
  geom_boxplot(width = 0.1, color = "grey", alpha=0.5)+
  theme(legend.position = "none",
        plot.title = element_text(size = 11))+
  ggtitle("GAPDH log-normalized counts")+
  geom_jitter(height = 0, width = 0.1,alpha=0.7)+
  scale_fill_viridis_d()+
  ylab("log10-normalized Counts") +
  xlab("Type")




ggscatter(imp, x="MelanA", y="GRSLND")

theme(axis.text.x=element_blank(),
      axis.ticks.x=c(0,50,100,150,200),
      axis.text.y = element_blank(),
      axis.ticks.y=c(0,5000,15000,25000))