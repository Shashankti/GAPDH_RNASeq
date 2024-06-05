#only TCGA LAML data
#overlapping with patient metadata

library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)


clinical <- GDCquery_clinic("TCGA-SKCM")
any(colnames(clinical) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical[,c(9,38,44)]

clinical$deceased <- ifelse(clinical$vital_status=="Alive",FALSE,TRUE)
clinical$overall_survival <- ifelse(clinical$vital_status=="Alive",
                                    clinical$days_to_last_follow_up,
                                    clinical$days_to_death)
clinical <- clinical[!is.na(clinical$overall_survival),]


query_TCGA <- GDCquery(project = "TCGA-LAML",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts",
                       experimental.strategy = "RNA-Seq")
GDCdownload(query_TCGA, method = "api", files.per.chunk = 100,
            directory = "~/TCGA_Data/TCGA-LAML/")
tcga_laml <- GDCprepare(query_TCGA, directory = "~/TCGA_Data/TCGA-LAML/",summarizedExperiment = TRUE)
#counts matrix
TCGA_df_2 <- assay(tcga_laml)
gene_metadata <- as.data.frame(rowData(tcga_laml))

#coldata
coldata_clinical <- as.data.frame(colData(tcga_laml))

dds_clinical <- DESeqDataSetFromMatrix(countData = TCGA_df_2,
                                       colData = coldata_clinical,
                                       design = ~1)

keep <- rowSums(counts(dds_clinical)) >= 10
dds_clinical <- dds_clinical[keep,]

#vst 
vsd_clinical <- vst(dds_clinical,blind = FALSE)
vsd_matrix <- assay(vsd_clinical)
#gapdh
laml_gapdh <- vsd_matrix %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id)%>%
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "GAPDH")

median_value <- median(laml_gapdh$counts)

laml_gapdh$strata <- ifelse(laml_gapdh$counts >= median_value,"HIGH","LOW")
laml_gapdh$case_id <- gsub('-03.*', '', laml_gapdh$case_id)
laml_gapdh <- merge(laml_gapdh, clinical, by.x = 'case_id', by.y = 'submitter_id')


# fitting survival curve -----------
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = laml_gapdh)
fit
surv_plots <- 
  ggsurvplot(fit,conf.int = FALSE,
           data = laml_gapdh,
           legend.labs=c("High","Low"),
           legend.title="GAPDH",
           title="Kaplan-Meier Curve for AML Survival",
           pval = T,
           risk.table = T,
           linetype = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           palette = c("#990000", "#000099"))
tiff("GAPDH_ANL_Survival_wo_error.tiff",width = 35,height = 21,units = 'cm',res = 300)
surv_plots
dev.off()

fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = laml_gapdh)
fit2
pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
print(pv)

##---------------------------------------------------------------##
# SKCM Survival

#load data from tcga-skcm
mrna_query <- GDCquery(project = "TCGA-SKCM",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts",
                       experimental.strategy = "RNA-Seq")


GDCdownload(mrna_query, method = "api", files.per.chunk = 100,
            directory = "~/TCGA_Data/TCGA-SKCM/")
mrna_df <- GDCprepare(mrna_query, directory = "~/TCGA_Data//TCGA-SKCM/")
# extract what you need from mrna_df$ to add to metadata
mrna_meta <- mrna_df$sample
mrna_meta <- cbind(mrna_meta, mrna_df$definition)
mrna_meta <- cbind(mrna_meta,  mrna_df$barcode)
skcm_df <- assay(mrna_df)

gene_metadata <- as.data.frame(rowData(mrna_df))

#coldata
coldata_clinical <- as.data.frame(colData(mrna_df))

dds_clinical <- DESeqDataSetFromMatrix(countData = skcm_df,
                                       colData = coldata_clinical,
                                       design = ~1)

keep <- rowSums(counts(dds_clinical)) >= 10
dds_clinical <- dds_clinical[keep,]

#vst 
vsd_clinical <- vst(dds_clinical,blind = FALSE)
vsd_matrix <- assay(vsd_clinical)

#gapdh
laml_gapdh <- vsd_matrix %>%
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id)%>%
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "GRASLND")

median_value <- median(laml_gapdh$counts)

laml_gapdh$strata <- ifelse(laml_gapdh$counts >= median_value,"HIGH","LOW")
laml_gapdh$case_id <- gsub('-06.*|-07.*|-01.*|-11.*', '', laml_gapdh$case_id)
laml_gapdh <- merge(laml_gapdh, clinical, by.x = 'case_id', by.y = 'submitter_id')
median_value <- median(laml_gapdh$counts)
laml_gapdh$strata <- ifelse(laml_gapdh$counts >= median_value,"HIGH","LOW")


# fitting survival curve -----------
fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = laml_gapdh)
fit
surv_plots <- 
  ggsurvplot(fit,conf.int = FALSE,
             data = laml_gapdh,
             legend.labs=c("High","Low"),
             legend.title="GRASLND",
             title="Kaplan-Meier Curve for SKCM Survival",
             pval = T,
             risk.table = T,
             linetype = "strata",
             surv.median.line = "hv",
             ggtheme = theme_bw(),
             palette = c("#990000", "#000099"))
tiff("../SKCM/Survival_plot_wo_error.tiff",width = 35,height = 21,units = 'cm',res = 300)
surv_plots
dev.off()

fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = laml_gapdh)
fit2
pv <- ifelse ( is.na(fit2),next,(round(1 - pchisq(fit2$chisq, length(fit2$n) - 1),3)))[[1]]
print(pv)
