#DESeq2 
#simple analysis with 1 factors and 2 levels
#load libraries
library (DESeq2)
library (tidyverse)
library (airway)

# Step 1: preparing count data ----------------
# read in counts data
counts_data <- read.table ('./Data/genes_numreads.tsv')
head(counts_data)

# read in sample info
colData<-read.table ('./Data/sample_info_all.tsv')
view (colData)

#removing columns in counts data to make specific count data for each condition
#
# count data for Norm shRNA vs control
counts_data_Norm = subset (counts_data, select = -c( Molm13.low.Rep1 , Molm13.low.Rep2 , Molm13.Gal.Rep1 , Molm13.Gal.Rep2 , Molm13.shRNA1.Low , Molm13.shRNA1.Gal ,Molm13.shRNA7.Low ,  Molm13.shRNA7.Gal ))
head (counts_data_Norm)
#
#
# count data for Low shRNA vs control
counts_data_Low = subset (counts_data, select = -c( Molm13.Norm.Rep1 , Molm13.Norm.Rep2 ,  Molm13.Gal.Rep1 , Molm13.Gal.Rep2 , Molm13.shRNA1.Norm , Molm13.shRNA1.Gal , Molm13.shRNA7.Norm , Molm13.shRNA7.Gal ))
head (counts_data_Low)
#
#
# count data for Gal shRNA vs control
counts_data_Gal = subset (counts_data, select = -c( Molm13.Norm.Rep1	, Molm13.Norm.Rep2	, Molm13.low.Rep1	, Molm13.low.Rep2	, Molm13.shRNA1.Norm	, Molm13.shRNA1.Low	,	Molm13.shRNA7.Norm , Molm13.shRNA7.Low ))
head (counts_data_Gal)
#
#
#count data for parental Low vs Norm 
counts_data_PLN = subset (counts_data, select = -c( Molm13.Gal.Rep1 , Molm13.Gal.Rep2 , Molm13.shRNA1.Norm	, Molm13.shRNA1.Low	, Molm13.shRNA1.Gal , Molm13.shRNA7.Norm , Molm13.shRNA7.Low , Molm13.shRNA7.Gal ))
head (counts_data_PLN)
#
#
#count data for parental Gal vs Norm
counts_data_PGN = subset (counts_data, select = -c (Molm13.low.Rep1	, Molm13.low.Rep2 ,  Molm13.shRNA1.Norm	, Molm13.shRNA1.Low	, Molm13.shRNA1.Gal , Molm13.shRNA7.Norm , Molm13.shRNA7.Low , Molm13.shRNA7.Gal))
head (counts_data_PGN )
#
#
#count data for parental Gal vs Low 
counts_data_PGL = subset (counts_data, select = -c ( Molm13.Norm.Rep1	, Molm13.Norm.Rep2,  Molm13.shRNA1.Norm	, Molm13.shRNA1.Low	, Molm13.shRNA1.Gal , Molm13.shRNA7.Norm , Molm13.shRNA7.Low , Molm13.shRNA7.Gal))
head (counts_data_PGL)
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
dds <- DESeqDataSetFromMatrix (countData = round(counts_data), colData = colData , design = ~ condition)
ddsNorm <- DESeqDataSetFromMatrix (countData = round(counts_data_Norm), colData = colDataNorm , design = ~ condition)
ddsLow <- DESeqDataSetFromMatrix (countData = round(counts_data_Low), colData = colDataLow , design = ~ condition)
ddsGal <- DESeqDataSetFromMatrix (countData = round(counts_data_Gal), colData = colDataGal , design = ~ condition)
ddsPLN <- DESeqDataSetFromMatrix (countData = round(counts_data_PLN), colData = colDataPLN , design = ~ condition)
ddsPGN <- DESeqDataSetFromMatrix (countData = round(counts_data_PGN), colData = colDataPGN , design = ~ condition)
ddsPGL <- DESeqDataSetFromMatrix (countData = round(counts_data_PGL), colData = colDataPGL , design = ~ condition)

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total

keep <- rowSums(counts(dds)) >= 10
dds <- dds [keep,]
dds

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
ddsNorm <- DESeq(ddsNorm)
ddsLow<- DESeq(ddsLow)
ddsGal <- DESeq(ddsGal)
ddsPLN <- DESeq(ddsPLN)
ddsPGN <- DESeq(ddsPGN)
ddsPGL <- DESeq(ddsPGL)
#
res_all <- results (dds)
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




