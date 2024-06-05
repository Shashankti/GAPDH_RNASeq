#run analysis with only selected samples
final_merged_counts
col_data_2 <- data.table::fread("Sample_info_2.csv",header = F)

col_data_2 <- col_data_2%>% remove_rownames() %>% column_to_rownames(var = "V1")
colnames(col_data_2) <- c("Source", "Type")
#col_data <- col_data[which(rownames(col_data) %in% colnames(final_merged_counts)),]

`%nin%` = Negate(`%in%`)

#removing extra samples
diff_counts <- final_merged_counts[,rownames(col_data_2)]

#making sure the row names in colData matches the column names in counts_data 
all(colnames(final_merged_counts) %in% rownames(col_data))

#check the order
all(colnames(final_merged_counts) == rownames(col_data))

# Step 2: construct a DESeqDataSet object
dds_2f <- DESeqDataSetFromMatrix (countData = diff_counts, colData = col_data_2 , design = ~ Type)

# keep_2f <- rowSums(counts(dds_2f)) >= 10
# dds_2f <- dds_2f [keep_2f,]
# dds_2f

#get normalised counts information
dds_2f <- estimateSizeFactors(dds_2f)
sizeFactors(dds_2f)
dds_2f$Type <- relevel(dds_2f$Type,ref="Control")

# Step 3: Run DESeq ----------------------

dds_2f <- DESeq(dds_2f)
resultnames_2f <- resultsNames(dds_2f)
resultnames_2f

result_TCGA <- results(dds_2f, contrast = c("Type", "Tumor", "Control")) 
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


## working with TCGA data
library(recount3)
library(interactiveDisplayBase)

human_projects <- available_projects()
subset(human_projects, file_source == "gtex" & project_type == "data_sources")
selected_study <- interactiveDisplayBase::display(human_projects)

proj_info <- subset(
  human_projects,
  project="BONE_MARROW" & project_type =="data_sources"
)
proj_info <- human_projects[human_projects$project=="BONE_MARROW",]
res_bone_marrow <- create_rse(proj_info)
bm_counts <- assay(res_bone_marrow)
gtcga <- na.omit(merge(TCGA_df, bm_counts, by="row.names", all = TRUE))
gtcga <- gtcga %>% remove_rownames() %>% column_to_rownames(var = "Row.names")