#running edgeR and Wilcoxons rank sum test
#edgeR TMM normalize
conditions <- as.data.frame(t(col_data))[2,]
conditions <- factor(t(conditions))
y <- DGEList(counts = na.omit(final_merged_counts[non_canonical_list$V2,]), group = conditions)
## Remove rows conssitently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep, keep.lib.sizes = FALSE]
## Perform TMM normalization and convert to CPM (Counts Per Million)
y <- calcNormFactors(y, method = "TMM")
count_norm <- cpm(y)
count_norm <- as.data.frame(count_norm)
count_norm$gene_name <- gene_IDs[match(rownames(count_norm), gene_IDs[,1]),2]
count_norm2 <- count_norm %>% remove_rownames %>%
  column_to_rownames(var = "gene_name")
count_norm2 <- count_norm2[rownames(count_norm2) %in% non_canonical_list$V1,]
# Run the Wilcoxon rank-sum test for each gene
count_norm$gene_name <- NULL
pvalues <- sapply(1:nrow(count_norm), function(i){
  data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions)
  p <- wilcox.test(gene~conditions, data)$p.value
  return(p)})
fdr <- p.adjust(pvalues, method = "fdr")
# Calculate the fold-change for each gene
conditionsLevel <- levels(conditions)
dataCon1 <- count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2 <- count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
# Output results based on the FDR threshold 0.05
outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
rownames(outRst) <- rownames(count_norm)
outRst <- na.omit(outRst)
fdrThres <- 0.05
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(outRst), mart= mart)
outRst$GeneName <- gene_IDs[match(rownames(outRst), gene_IDs[,1]),2]


EnhancedVolcano::EnhancedVolcano(outRst[outRst$GeneName %in% non_canonical_list$V1,],
                                 lab = outRst[outRst$GeneName %in% non_canonical_list$V1,]$GeneName,
                                          x = 'log2foldChange',
                                          y = 'FDR',
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

##make violin plot


##-------------------------------------------------#
#Looking at possible batch correction methods

library(sva)