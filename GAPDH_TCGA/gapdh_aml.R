#add GTex data to final merged counts

#read in the gtex data and cleanup
normal_blood <- data.table::fread("~/Downloads/gene_reads_2017-06-05_v8_spleen.gct.gz")
normal_blood$Name <- sub('\\.[0-9]*$','',normal_blood$Name)
normal_blood$id <- NULL
normal_blood$Description <- NULL
normal_blood <- normal_blood %>% remove_rownames() %>% column_to_rownames(var = "Name")

#remove normal samples from tcga skcm data
`%nin%` = Negate(`%in%`)
# mrna_meta <- as.data.table(mrna_meta)
# mrna_meta_2 <- mrna_meta[V2 %nin% c("Additional Metastatic","Solid Tissue Normal")]
# mrna_counts <- mrna_df[,colnames(mrna_df) %in% mrna_meta_2$V3]
rownames(mrna_counts) <- sub('\\.[0-9]*$','',rownames(mrna_counts))

gtex_merged_counts <- na.omit(merge(final_merged_counts, normal_blood, by="row.names", all = TRUE))
gtex_merged_counts <- gtex_merged_counts %>% remove_rownames() %>% column_to_rownames(var = "Row.names")
#final_merged_counts <- final_merged_counts[,which(colnames(final_merged_counts) %in% colnames(file2))]
gtex_coldata <- data.frame(ID=colnames(normal_blood),
                           Source = "GTEX",
                           Type="Control")
gtex_coldata <- gtex_coldata %>% remove_rownames() %>% column_to_rownames(var = "ID")
gtex_coldata <- rbind(col_data,gtex_coldata)

conditions_aml <- factor(as.data.frame(gtex_coldata)[,2])
y <- DGEList(counts =gtex_merged_counts, group = conditions_aml)
## Remove rows conssitently have zero or very low counts
keep <- filterByExpr(y)
y <- y[keep, keep.lib.sizes = FALSE]
## Perform TMM normalization and convert to CPM (Counts Per Million)
y <- calcNormFactors(y, method = "TMM")
count_norm <- cpm(y)
count_norm <- as.data.frame(count_norm)
rownames(count_norm) <- gsub("\\..*","",rownames(count_norm))
count_norm$gene_name <- gene_IDs[match(rownames(count_norm), gene_IDs[,1]),2]
count_norm2 <- count_norm %>% remove_rownames %>%
  column_to_rownames(var = "gene_name")
count_norm2 <- count_norm2[rownames(count_norm2) %in% non_canonical_list$V1,]
# Run the Wilcoxon rank-sum test for each gene
count_norm$gene_name <- NULL
pvalues <- sapply(1:nrow(count_norm), function(i){
  data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), conditions_aml)
  p <- wilcox.test(gene~conditions_aml, data)$p.value
  return(p)})
fdr <- p.adjust(pvalues, method = "fdr")



# Calculate the fold-change for each gene
conditionsLevel <- levels(conditions_aml)
dataCon1 <- count_norm[,c(which(conditions_aml==conditionsLevel[1]))]
dataCon2 <- count_norm[,c(which(conditions_aml==conditionsLevel[2]))]
foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))
# Output results based on the FDR threshold 0.05
outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)
rownames(outRst) <- rownames(count_norm)
outRst <- na.omit(outRst)
fdrThres <- 0.05
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(count_norm), mart= mart)
outRst$GeneName <- gene_IDs[match(rownames(outRst), gene_IDs[,1]),2]

gapdh <- count_norm["ENSG00000111640",]
gapdh <- as.data.frame(t(gapdh))
gapdh$type <- y$samples[match(rownames(gapdh), rownames(y$samples)),1]
colnames(gapdh)[1] <- "Count"


viol_plot <- ggplot(gapdh, aes(x = type, y = log(Count+1,base = 2), fill = type))+
  geom_violin(trim = FALSE, width=0.5,)+
  theme_cowplot()+
  geom_boxplot(width = 0.15, color = "black", alpha=0.5)+
  theme(legend.position = "none",
        plot.title = element_text(size = 11))+
  ggtitle("GRASLND Normalized counts")+
  geom_jitter(width = 0.05,alpha=0.4,size=0.2)+
  scale_fill_manual(values = c("#a0044d","#610051"))+
  ylab("Expression - log2(TPM+1)") +
  xlab("SKCM
       (num(N)=701;num(T)=471)")+
  geom_signif(comparisons = list(c("Normal","Tumor")),
              map_signif_level = T,
              test = wilcox.test,
              y_position = 6.8, tip_length = 0.1, vjust = 0.2)
  annotate("text",x = 0.58, xend = 1.9, y= 7, label = "Mean LFC = 4.56")+
  annotate("text",x = 0.58, xend = 1.9, y= 6.5, label = "pvalue = 3e-153")

png("~/Documents/SKCM/viol_plot.png",width = 35,height = 21,units = 'cm',res = 300)
plot(viol_plot)
dev.off()

stat.test <- gapdh %>%
  (Count ~ type) %>%
  add_significance()
stat.test
fold.change <- log(mean(grslnd[grslnd$type=="Tumor",1])/mean(grslnd[grslnd$type=="Normal",1]),base=2)
fold.change <- log(median(grslnd[grslnd$type=="Tumor",1])/median(grslnd[grslnd$type=="Normal",1]),base=2)





