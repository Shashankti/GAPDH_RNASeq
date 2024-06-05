#DESeq2 
#simple analysis with one condition and 2 factors
#load libraries
library (DESeq2)
library (tidyverse)
library (airway)

# Step 1: preparing count data ----------------
# read in counts data
counts_data <- read.table ('./Data/genes_tpm.tsv',header = TRUE)
head(counts_data)
counts_data <- setNames(counts_data,c("Molm13.Norm.Rep1","Molm13.Norm.Rep2","Molm13.low.Rep1","Molm13.low.Rep2","Molm13.Gal.Rep1","Molm13.Gal.Rep2","Molm13.shRNA1.Norm","Molm13.shRNA1.Low","Molm13.shRNA1.Gal","Molm13.shRNA7.Norm","Molm13.shRNA7.Low","Molm13.shRNA7.Gal"))

# read in sample info
colData<-read.table ('./Data/sample_info_all.tsv')
colData_grouped <- read.table('./Data/sample_info_all_grouped.tsv')
view (colData)

#removing columns in counts data to make specific count data for each condition
#
# count data for Norm shRNA vs control
counts_data_Norm = counts_data[,c(1,2,7,10)]
#
#
# count data for Low shRNA vs control
counts_data_Low = counts_data[,c(3,4,8,11)]
#
# count data for Gal shRNA vs control
counts_data_Gal = counts_data[,c(5,6,9,12)]
#
#
#count data for parental Low vs Norm 
counts_data_PLN = counts_data[,c(1,2,3,4)]
#
#
#count data for parental Gal vs Norm
counts_data_PGN = counts_data[,c(1,2,5,6)]
#
#
#count data for parental Gal vs Low 
counts_data_PGL = counts_data[,c(3,4,5,6)]
#
#
#read new sample info that matches our count data 

colDataNorm <- read.table ('./Data/sample_info_norm.tsv')
colDataLow  <- read.table ('./Data/sample_info_low.tsv')
colDataGal <- read.table ('./Data/sample_info_gal.tsv')
colDataPLN <- read.table ('./Data/sample_info_PLN.tsv')
colDataPGN <- read.table ('./Data/sample_info_PGN.tsv')
colDataPGL <- read.table ('./Data/sample_info_PGL.tsv')

#making sure the row names in colData matches the column names in counts_data 

all(colnames(counts_data) %in% rownames(colData))
all(colnames(counts_data_Norm) %in% rownames(colDataNorm))
all(colnames(counts_data_Low) %in% rownames(colDataLow))
all(colnames(counts_data_Gal) %in% rownames(colDataGal))
all(colnames(counts_data_PLN) %in% rownames(colDataPLN))
all(colnames(counts_data_PGN) %in% rownames(colDataPGN))
all(colnames(counts_data_PGL) %in% rownames(colDataPGL))

#check the order
all(colnames(counts_data) == rownames(colData))
all(colnames(counts_data_Norm) == rownames(colDataNorm))
all(colnames(counts_data_Low) == rownames(colDataLow))
all(colnames(counts_data_Gal) == rownames(colDataGal))
all(colnames(counts_data_PLN) == rownames(colDataPLN))
all(colnames(counts_data_PGN) == rownames(colDataPGN))
all(colnames(counts_data_PGL) == rownames(colDataPGL))

# Step 2: construct a DESeqDataSet object
dds_copy <- DESeqDataSetFromMatrix (countData = round(counts_data), colData = colData_grouped , design = ~ condition)
ddsNorm <- DESeqDataSetFromMatrix (countData = round(counts_data_Norm), colData = colDataNorm , design = ~ condition)
ddsLow <- DESeqDataSetFromMatrix (countData = round(counts_data_Low), colData = colDataLow , design = ~ condition)
ddsGal <- DESeqDataSetFromMatrix (countData = round(counts_data_Gal), colData = colDataGal , design = ~ condition)
ddsPLN <- DESeqDataSetFromMatrix (countData = round(counts_data_PLN), colData = colDataPLN , design = ~ condition)
ddsPGN <- DESeqDataSetFromMatrix (countData = round(counts_data_PGN), colData = colDataPGN , design = ~ condition)
ddsPGL <- DESeqDataSetFromMatrix (countData = round(counts_data_PGL), colData = colDataPGL , design = ~ condition)

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total

keep <- rowSums(counts(dds_copy)) >= 10
dds_copy <- dds_copy [keep,]
dds_copy

##
keepNorm <- rowSums(counts(ddsNorm)) >= 10
ddsNorm <- ddsNorm [keepNorm]
ddsNorm
##
keepLow <- rowSums(counts(ddsLow)) >= 10
ddsLow <- ddsLow [keepLow]
ddsLow
##
keepGal <- rowSums(counts(ddsGal)) >= 10
ddsGal <- ddsGal [keepGal]
ddsGal
##
keepPLN<- rowSums(counts(ddsPLN)) >= 10
ddsPLN <- ddsPLN[keepPLN]
ddsPLN
##
keepPGN <- rowSums(counts(ddsPGN)) >= 10
ddsPGN <- ddsPGN[keepPGN]
ddsPGN
##
keepPGL <- rowSums(counts(ddsPGL)) >= 10
ddsPGL <- ddsPGL[keepPGL]
ddsPGL

# set the factor level only for parental samples 
ddsPLN$condition <- relevel(ddsPLN$condition, ref = "norm")
ddsPGN$condition <- relevel(ddsPGN$condition, ref = "norm")
ddsPGL$condition <- relevel(ddsPGL$condition, ref = "low")

# Step 3: Run DESeq ----------------------

dds <- DESeq(dds)
dds_copy <- DESeq(dds_copy)
ddsNorm <- DESeq(ddsNorm)
ddsLow<- DESeq(ddsLow)
ddsGal <- DESeq(ddsGal)
ddsPLN <- DESeq(ddsPLN)
ddsPGN <- DESeq(ddsPGN)
ddsPGL <- DESeq(ddsPGL)
#
res_all <- results (dds_copy)
res_all
#
res_shNorm_vs_ctNorm <- results(ddsNorm)
res_shNorm_vs_ctNorm  
#
res_shLow_vs_ctLow <- results(ddsLow)
res_shLow_vs_ctLow
#
res_shGal_vs_ctGal <- results(ddsGal)
res_shGal_vs_ctGal
#
res_pLow_vs_pNorm <- results (ddsPLN)
res_pLow_vs_pNorm
# 
res_pGal_vs_pNorm <- results (ddsPGN)
res_pGal_vs_pNorm
# 
res_pGal_vs_pLow <- results (ddsPGL)
res_pGal_vs_pNorm

# Explore Results ----------------
summary (res_all)
summary (res_shNorm_vs_ctNorm)
summary (res_shLow_vs_ctLow)
summary (res_shGal_vs_ctGal)
summary (res_pLow_vs_pNorm)
summary (res_pGal_vs_pNorm)
summary (res_pGal_vs_pLow)

# exporting results to CSV files
write.csv(as.data.frame (res_all) , file = "./Analysis/DEseq2_shRNA_vs_control_filtered.csv")
write.csv(as.data.frame (res_shNorm_vs_ctNorm) , file = "./Analysis/DEseq2_shNorm_vs_ctNorm_filtered.csv")
write.csv(as.data.frame (res_shLow_vs_ctLow) , file = "./Analysis/DEseq2_shLow_vs_ctLow_filtered.csv")
write.csv(as.data.frame (res_shGal_vs_ctGal) , file = "./Analysis/DEseq2_shGal_vs_ctGal_filtered.csv")
write.csv(as.data.frame (res_pLow_vs_pNorm) , file = "./Analysis/DEseq2_pLow_vs_pNorm_filtered.csv")
write.csv(as.data.frame (res_pGal_vs_pNorm) , file = "./Analysis/DEseq2_pGal_vs_pNorm_filtered.csv")
write.csv(as.data.frame (res_pGal_vs_pLow) , file = "./Analysis/DEseq2_pGal_vs_pLow_filtered.csv")


design(dds_copy) <- formula(~ type+condition)
dds_copy <- DESeq(dds_copy)
res_mf <- results(dds_copy,
                  contrast = c("type","Gal","Low"))
head(res_mf)

dds_copy <- estimateSizeFactors(dds_copy)
normalised_counts <- counts(dds_copy,normalized=TRUE)
write.table(normalised_counts,file = "./Analysis/All_normalised_counts.txt",sep = "\t",quote =F,col.names = NA)

###
ddsGal <- estimateSizeFactors(ddsGal)
normalised_counts_Gal <- counts(ddsGal,normalized=TRUE)
write.table(normalised_counts_Gal,file = "./Analysis/Gal_normalised_counts.txt",sep="\t",quote=F,col.names = NA)
##
ddsNorm <- estimateSizeFactors(ddsNorm)
normalised_counts_Norm <- counts(ddsNorm,normalized=TRUE)
write.table(normalised_counts_Norm,file = "./Analysis/Norm_normalised_counts.txt",sep="\t",quote=F,col.names = NA)
##
ddsLow <- estimateSizeFactors(ddsLow)
normalised_counts_low <- counts(ddsLow,normalized=TRUE)
write.table(normalised_counts_low,file = "./Analysis/low_normalised_counts.txt",sep="\t",quote=F,col.names = NA)
##
ddsPLN <- estimateSizeFactors(ddsPLN)
normalised_counts_PLN <- counts(ddsPLN,normalized=TRUE)
write.table(normalised_counts_PLN,file = "./Analysis/PLN_normalised_counts.txt",sep="\t",quote=F,col.names = NA)
##
ddsPGN <- estimateSizeFactors(ddsPGN)
normalised_counts_PGN <- counts(ddsPGN,normalized=TRUE)
write.table(normalised_counts_PGN,file = "./Analysis/PGN_normalised_counts.txt",sep="\t",quote=F,col.names = NA)
##
ddsPGL <- estimateSizeFactors(ddsPGL)
normalised_counts_PGL <- counts(ddsPGL,normalized=TRUE)
write.table(normalised_counts_PGL,file = "./Analysis/PGL_normalised_counts.txt",sep="\t",quote=F,col.names = NA)

