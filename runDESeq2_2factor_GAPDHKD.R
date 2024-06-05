#DESeq2 
#analysis with complex contrast (2 factor with interactions)
# tutorial: https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
#load libraries
library (DESeq2)
library (tidyverse)
library (airway)

# Step 1: preparing count data ----------------
# read in counts data
counts_data <- read.table ('./Data/genes_numreads.tsv')
head(counts_data)

# read in sample info
colData_2f<-read.table ('./Data/sample_info.tsv')
#view (colData_2f)



#making sure the row names in colData matches the column names in counts_data 

all(colnames(counts_data) %in% rownames(colData_2f))

#check the order
all(colnames(counts_data) == rownames(colData_2f))


#Making contrasts
#1)get the model matrix
mod_mat_2f <- model.matrix(design(dds_2f), colData(dds_2f))

# Step 2: construct a DESeqDataSet object
dds_2f <- DESeqDataSetFromMatrix (countData = round(counts_data), colData = colData_2f , design = ~ condition + type + type:condition )
dds_1f <- DESeqDataSetFromMatrix (countData = round(counts_data[,1:6]), colData = colData_2f[1:6,] , design = ~ condition)

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total

keep_2f <- rowSums(counts(dds_1f)) >= 10
dds_1f <- dds_1f [keep_2f,]
dds_1f

#get normalised counts information
dds_2f <- estimateSizeFactors(dds_2f)
sizeFactors(dds_2f)



#releveal condition
dds_2f$condition <- relevel(dds_2f$condition, ref = "Norm")
dds_1f$condition <- relevel(dds_1f$condition, ref = "Norm")


# Step 3: Run DESeq ----------------------

dds_2f <- DESeq(dds_2f)
resultnames_2f <- resultsNames(dds_2f)
resultnames_2f 

dds_1f <- DESeq(dds_1f)
resultnames_2f <- resultsNames(dds_1f)




# get the model matrix
mod_mat_2f <- model.matrix(design(dds_2f), colData(dds_2f))
mod_mat_2f

# Define coefficient vectors for each condition
norm_parental <- colMeans(mod_mat_2f[dds_2f$condition == "Norm" & dds_2f$type == "parental", ])
low_parental <- colMeans(mod_mat_2f[dds_2f$condition == "Low" & dds_2f$type == "parental", ])
gal_parental <- colMeans(mod_mat_2f[dds_2f$condition == "Gal" & dds_2f$type == "parental", ])
norm_shRNA <- colMeans(mod_mat_2f[dds_2f$condition == "Norm" & dds_2f$type == "shRNA", ])
low_shRNA <- colMeans(mod_mat_2f[dds_2f$condition == "Low" & dds_2f$type == "shRNA", ])
gal_shRNA <- colMeans(mod_mat_2f[dds_2f$condition == "Gal" & dds_2f$type == "shRNA", ])

#define contrast of interest from  vectors
res_shNorm_vs_ctNorm_2 <- results(dds_2f,contrast = list(c("Intercept","type_shRNA_vs_parental") ))

res_pLow_vs_pNorm_2 <-results (dds_2f, contrast = low_parental - norm_parental) 
res_pGal_vs_pNorm_2 <- results(dds_2f, contrast = gal_parental - norm_parental)
res_pGal_vs_pLow_2 <- results(dds_2f, contrast = gal_parental - low_parental)
res_shNorm_vs_ctNorm_2 <- results(dds_2f, contrast = norm_shRNA - norm_parental)
res_shLow_vs_ctLow_2 <- results(dds_2f, contrast = low_shRNA - low_parental)
res_shGal_vs_ctGal_2 <- results(dds_2f, contrast = gal_shRNA - gal_parental)


res_control_vs_sh <- results(dds_2f,contrast = c('type','shRNA','parental'))


res_shNorm_vs_ctNorm_Unshrunk <- results(dds_2f, contrast = norm_shRNA - norm_parental,independentFiltering = T)
deRes <- as.data.frame(res_shNorm_vs_ctNorm_Unshrunk)
deRes <- na.omit(deRes)

deRes$baseMeanGroup <- groups_by_filter(deRes$baseMean, 8)

ggplot(deRes, aes(x=pvalue)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ baseMeanGroup, nrow = 2)

# Explore Results ----------------
summary (res_shNorm_vs_ctNorm_2)
summary (res_shLow_vs_ctLow_2)
summary (res_shGal_vs_ctGal_2)
summary (res_pLow_vs_pNorm)
summary (res_pGal_vs_pNorm)
summary (res_pGal_vs_pLow)


# Annotate Gene names
res_pLow_vs_pNorm_2 <-arrange(as.data.frame(res_pLow_vs_pNorm_2),padj) %>% rownames_to_column("GeneID")
res_pGal_vs_pNorm_2 <- arrange(as.data.frame(res_pGal_vs_pNorm_2),padj) %>% rownames_to_column("GeneID")
res_pGal_vs_pLow_2 <- arrange(as.data.frame(res_pGal_vs_pLow_2),padj) %>% rownames_to_column("GeneID")
res_shNorm_vs_ctNorm_2 <- arrange(as.data.frame(res_shNorm_vs_ctNorm_2),padj) %>% rownames_to_column("GeneID")
res_shLow_vs_ctLow_2 <- arrange(as.data.frame(res_shLow_vs_ctLow_2),padj) %>% rownames_to_column("GeneID")
res_shGal_vs_ctGal_2 <- arrange(as.data.frame(res_shGal_vs_ctGal_2),padj) %>% rownames_to_column("GeneID")

#pLow vs pNorm
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = res_pLow_vs_pNorm_2[,1], mart= mart)

res_pLow_vs_pNorm_2 <- arrange(merge(res_pLow_vs_pNorm_2,gene_IDs,by.x="GeneID",by.y="ensembl_gene_id"),padj)

# pGal_pNorm
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = res_pGal_vs_pNorm_2[,1], mart= mart)
res_pGal_vs_pNorm_2 <- arrange(merge(res_pGal_vs_pNorm_2,gene_IDs,by.x="GeneID",by.y="ensembl_gene_id"),padj)

#pGal vs pLow
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = res_pGal_vs_pLow_2[,1], mart= mart)

res_pGal_vs_pLow_2 <- arrange(merge(res_pGal_vs_pLow_2,gene_IDs,by.x="GeneID",by.y="ensembl_gene_id"),padj)

#shNorm vs ctNorm
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = res_shNorm_vs_ctNorm_2[,1], mart= mart)

res_shNorm_vs_ctNorm_2 <- arrange(merge(res_shNorm_vs_ctNorm_2,gene_IDs,by.x="GeneID",by.y="ensembl_gene_id"),padj)

# shLow vs ctLow

gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = res_shLow_vs_ctLow_2[,1], mart= mart)

res_shLow_vs_ctLow_2 <- arrange(merge(res_shLow_vs_ctLow_2,gene_IDs,by.x="GeneID",by.y="ensembl_gene_id"),padj)

#shGal vs ctGAl
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = res_shGal_vs_ctGal_2[,1], mart= mart)

res_shGal_vs_ctGal_2 <- arrange(merge(res_shGal_vs_ctGal_2,gene_IDs,by.x="GeneID",by.y="ensembl_gene_id"),padj)


#Get normalized counts
normalized_counts <- counts(dds_2f,normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts) %>% rownames_to_column("GeneID")
write.table(normalized_counts,file = "./Analysis/Normalized_counts_2f.txt",sep = "\t",quote = F,col.names = NA)

# Add gene names to normalised counts
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = normalized_counts[,1], mart= mart)

#!normalized_counts <- merge(normalized_counts,gene_IDs, by.x="GeneID",by.y="ensembl_gene_id")
normalized_counts$GeneName <- gene_IDs[match(normalized_counts[,1],gene_IDs[,1]),2]

# Creating a detailed results table
res_pLow_vs_pNorm_OE <- lfcShrink(dds_2f, contrast = low_parental - norm_parental,res = res_pLow_vs_pNorm_2)
res_pGal_vs_pNorm_OE <- results(dds_2f, contrast = gal_parental - norm_parental)
res_pGal_vs_pLow_OE <- results(dds_2f, contrast = gal_parental - low_parental)
res_shNorm_vs_ctNorm_OE <- results(dds_2f, contrast = norm_shRNA - norm_parental)
res_shLow_vs_ctLow_OE <- results(dds_2f, contrast = low_shRNA - low_parental)
res_shGal_vs_ctGal_OE <- results(dds_2f, contrast = gal_shRNA - gal_parental)


#Plot the dispersion estimates
plotDispEsts(dds_2f)
plotMA()
top_genes <- res_shNorm_vs_ctNorm_2$GeneID[1:10]
vsd <- vst(dds_2f)
pheatmap::pheatmap(assay(vsd)[top_genes,],annotation_col = sampleInfo)
sampleInfo <- as.data.frame(colData(dds_2f)[,c("condition","type")])





