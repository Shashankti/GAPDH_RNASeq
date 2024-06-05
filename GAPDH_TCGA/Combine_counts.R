#make new TCGA

#merge custom counts file into one matrix
library(plyr)
library(edgeR)
library(Rsubread)
library(dplyr)
library(tidyverse)
library(biomaRt)

?featureCounts
options(timeout = '40000')


#read in the extra counts files

ex1 <- read.table("00AACC/GSE74246_RNAseq_All_Counts.txt", row.names = 1)
counts_list <- list.files("GSM/GSE63569_RAW/",pattern = ".txt",full.names = TRUE)
my_files <- lapply(counts_list, data.table::fread)
colnames = c("ENSEMBL ID")
for(i in seq_along(my_files)){
  colnames(my_files[[i]])[1] <- colnames
}


my_filenames <- c(my_files[[1]],my_files[[2]])
y <- readDGE(files = my_files)

my_files %>% reduce(left_join, by = "Counts")

GSE_df <- na.omit(Reduce(function(...) merge(..., by = "ENSEMBL ID", all=TRUE), my_files))
GSE_df <- GSE_df %>% remove_rownames() %>% column_to_rownames(var = "ENSEMBL ID")
GSE_df$ENSEMBL_ID <- rownames(GSE_df) 
#convert gene name to ensembl ID
## control vs d3
mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
gene_IDs <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(ex1), mart= mart)

ex1$ENSEMBL_ID <- gene_IDs[match(rownames(ex1),gene_IDs[,2]),1]
ex1 <- na.omit(ex1)
ex1 <- ex1 %>% remove_rownames() %>% column_to_rownames(var = "ENSEMBL_ID")

merged_counts <- merge(GSE_df, ex1, by = "ENSEMBL_ID", all = TRUE)
merged_counts <- na.omit(merged_counts)
write.table(merged_counts,file = "Merged_counts_AML_Extra.txt", quote = F, sep = "\t")



#load in the merged counts and the TCGA data together

merged_counts <- read.table("Merged_counts_AML_Extra.txt")
merged_counts <- merged_counts[!duplicated(merged_counts$ENSEMBL_ID),]
merged_counts <- merged_counts %>% remove_rownames() %>% column_to_rownames(var = "ENSEMBL_ID")

TCGA_df <- assay(data)

rownames(TCGA_df) <- sub('\\.[0-9]*$', '', rownames(TCGA_df))

final_merged_counts <- na.omit(merge(TCGA_df, merged_counts, by="row.names", all = TRUE))
final_merged_counts <- final_merged_counts %>% remove_rownames() %>% column_to_rownames(var = "Row.names")
#final_merged_counts <- final_merged_counts[,which(colnames(final_merged_counts) %in% colnames(file2))]
final_merged_counts <- final_merged_counts[,-(192:217)]
#read samples info

cd <- data.table::fread("Sample_info.tsv", header = TRUE)
ID = NULL
ID <- as.data.frame(colnames(final_merged_counts))
write.table(ID,file = "Sample_info.csv", quote = F, col.names = F, row.names = F)


