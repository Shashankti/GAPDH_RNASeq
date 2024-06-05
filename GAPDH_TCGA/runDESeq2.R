#running DESEQ2
library(DESeq2)
library(biomaRt)
col_data <- read.table("Sample_info.csv", sep = ",")
col_data <- col_data %>% remove_rownames() %>% column_to_rownames(var = "V1")
colnames(col_data) <- c("Source", "Type")
#col_data <- col_data[which(rownames(col_data) %in% colnames(final_merged_counts)),]

`%nin%` = Negate(`%in%`)


#making sure the row names in colData matches the column names in counts_data 

all(colnames(final_merged_counts) %in% rownames(col_data))

#check the order
all(colnames(final_merged_counts) == rownames(col_data))

# Step 2: construct a DESeqDataSet object
dds_2f <- DESeqDataSetFromMatrix (countData = round(final_merged_counts), colData = col_data , design = ~ Type + Source )

# keep_2f <- rowSums(counts(dds_2f)) >= 10
# dds_2f <- dds_2f [keep_2f,]
# dds_2f

#get normalised counts information
dds_2f <- estimateSizeFactors(dds_2f)
sizeFactors(dds_2f)
dds_2f$Type <- relevel(dds_2f$Type,ref="Normal")

# Step 3: Run DESeq ----------------------

dds_2f <- DESeq(dds_2f)
resultnames_2f <- resultsNames(dds_2f)
resultnames_2f

result_TCGA <- results(dds_2f, contrast = c("Type", "AML", "Normal")) 
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(result_TCGA), mart= mart)
result_TCGA$GeneName <- gene_IDs[match(rownames(result_TCGA), gene_IDs[,1]),2]
result_TCGA <- as.data.frame(result_TCGA)
result_TCGA$log2FoldChange <- as.numeric(result_TCGA$log2FoldChange)

write.table(result_TCGA, file = "All_DEG.tsv", quote = F, sep = "\t", row.names = T)
#read in non canonical genes
non_canonical_list <- read.table("~/Documents/AML/NonCanonical.txt", header = F)

#plots
plot1 <- EnhancedVolcano::EnhancedVolcano(result_TCGA[result_TCGA$GeneName %in% non_canonical_list$V1,],
                                 lab = result_TCGA[result_TCGA$GeneName %in% non_canonical_list$V1,]$GeneName,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 selectLab = c("GAPDH", "RPL18A","RPS10","RPL15","MRPL32"),
                                 title = 'Subsetting NonCanonical RBPs',
                                 subtitle = 'Differential Expressions',
                                 caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                                 pCutoff = 0.05,
                                 FCcutoff = 0.5,
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
                                 colConnectors = 'black')
tiff("NoncanonicalRBPs_Subset.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(plot1)
dev.off()
png("NoncanonicalRBPs_Subset.png",width = 35,height = 21,units = 'cm',res = 300)
plot(plot1)
dev.off()


#only non canonical DESEq2
counts_data <- final_merged_counts
counts_data$GeneName <- gene_IDs[match(rownames(counts_data), gene_IDs[,1]),2]
counts_data <- counts_data[counts_data$GeneName %in% non_canonical_list$V1,]
counts_data$GeneName <- NULL

dds_non_can <- DESeqDataSetFromMatrix(countData = round(counts_data), colData = col_data, design = ~ Type )

# keep_2f <- rowSums(counts(dds_non_can)) >= 10
# dds_non_can <- dds_non_can[keep_2f,]
# dds_non_can <- estimateSizeFactors(dds_non_can)
dds_non_can$Type <- relevel(dds_non_can$Type, ref = "Normal")

dds_non_can <- DESeq(dds_non_can)
result_non_can <- results(dds_non_can, contrast = c("Type","AML","Normal"))
result_non_can$GeneName <- gene_IDs[match(rownames(result_non_can), gene_IDs[,1]),2]
result_non_can_df <- as.data.frame(result_non_can)
write.table(result_non_can_df, file = "Only_ncRBPs_DEGresult.tsv", quote = F, row.names = T, sep = "\t")


plot2 <- EnhancedVolcano::EnhancedVolcano(result_non_can_df,
                                 lab = result_non_can_df$GeneName,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 selectLab = c("GAPDH", "RPL18A","RPS10","RPL15","MRPL32"),
                                 title = 'Only NonCanonical RBPs',
                                 subtitle = 'Differential Expressions',
                                 caption = bquote(~Log[2]~ 'foldchange cutoff 1; p-value cutoff,0.05'),
                                 pCutoff = 0.05,
                                 FCcutoff = 0.5,
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
                                 colConnectors = 'black')

tiff("NoncanonicalRBPs_complete.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(plot2)
dev.off()
png("NoncanonicalRBPs_complete.png",width = 35,height = 21,units = 'cm',res = 300)
plot(plot2)
dev.off()



## Function to split results
get_upregulated <- function(df){
  key <- intersect(rownames(df)[which(df$log2FoldChange>=1)],
                   rownames(df)[which(df$pvalue<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}

get_downregulated <- function(df){
  key <- intersect(rownames(df)[which(df$log2FoldChange<=-1)],
                   rownames(df)[which(df$pvalue<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}

up_reg <- get_upregulated(result_TCGA)


#take all GAPDH in normal vs cancer cells
GAPDH_all <- normalized_counts[rownames(normalized_counts)=="ENSG00000111640",]

normalized_counts <- counts(dds_2f,normalized=TRUE)
normalized_counts <- as.data.frame(normalized_counts) %>% rownames_to_column("GeneID")

d_mat <- as.matrix(t())

ggplot(GAPDH_all, aes())


##check for batch effects
set <- newSeqExpressionSet(counts(dds_2f))
idx  <- rowSums(counts(set) > 5) >= 2
set  <- set[idx, ]
set <- betweenLaneNormalization(set, which="upper")
not.sig <- rownames(result_TCGA)[which(result_TCGA$padj > .05)]
empirical <- rownames(set)[ rownames(set) %in% not.sig ]
set <- RUVg(set, empirical, k=2)
pData(set)
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(pData(set)[, i] ~ dds_2f$Source, vertical = TRUE, main = paste0("W", i))
  abline(h = 0)
}


ddsruv <- dds_2f
ddsruv$W1 <- set$W_1
ddsruv$W2 <- set$W_2
design(ddsruv) <- ~ W1 + W2 + Type + Source

# ddsruv <- DESeqDataSetFromMatrix(countData = counts(set), colData = pData(set), design = ~ W_1 + W_2+ x)

# keep_2f <- rowSums(counts(dds_non_can)) >= 10
# dds_non_can <- dds_non_can[keep_2f,]
# dds_non_can <- estimateSizeFactors(dds_non_can)
ddsruv$Type <- relevel(ddsruv$Type, ref = "Normal")

ddsruv <- DESeq(ddsruv)
result_ruv <- results(ddsruv, contrast = c("Type","AML","Normal"))
result_ruv$GeneName <- gene_IDs[match(rownames(result_ruv), gene_IDs[,1]),2]
result_ruv <- as.data.frame(result_ruv)

EnhancedVolcano::EnhancedVolcano(result_ruv[result_ruv$GeneName %in% non_canonical_list$V1,],
                                 lab = result_ruv[result_ruv$GeneName %in% non_canonical_list$V1,]$GeneName,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 selectLab = c("GAPDH", "RPL18A","RPS10","RPL15","MRPL32"),
                                 title = 'Subsetting NonCanonical RBPs',
                                 subtitle = 'Differential Expression',
                                 caption = bquote(~Log[2]~ 'foldchange cutoff 0.5; p-value cutoff,0.05'),
                                 pCutoff = 0.05,
                                 FCcutoff = 0.5,
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
                                 colConnectors = 'black')


