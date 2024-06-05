# Perform GO analysis on significantly expressed genes

library(clusterProfiler)
library(goseq)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggridges)
library(pathview)
library(ggplot2)
set.seed(2022)

#Select up regulated gene
up_genes <- rownames(sigOE_qValue[sigOE_qValue$log2FoldChange>0,])
#make the over representation test
GO_res <- enrichGO(gene = up_genes,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENSEMBL",
                   ont = "BP")
as.data.frame(GO_res)
# plot the results
fit1 <- plot(barplot(GO_res,showCategory = 30,font.size = 7,title = "Significantly up-regulated genes in Normal Glucose"))
fit1

# same process for the downregulated genes
down_genes <- rownames(sigOE_qValue[sigOE_qValue$log2FoldChange<0,])

GO_res_d <- enrichGO(gene = down_genes,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENSEMBL",
                   ont = "BP")
as.data.frame(GO_res_d)
fit2 <- plot(barplot(GO_res_d,showCategory = 30,font.size = 7,title = "Significantly down-regulated genes in Normal Glucose"))
fit2

GO_res_all <- enrichGO(gene=rownames(sigOE_qValue),
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP")

x1 <- pairwise_termsim(GO_res_d)
emapplot(x1,showCategory = 15)

# start the GSEA with DESEQ2 output
deRes_GSEA <- df_res_control_vs_shRNA[df_res_control_vs_shRNA$baseMean>50,]

#deRes_GSEA <- deRes[deRes$baseMean>50,]
deRes_GSEA <- deRes_GSEA[order(-deRes_GSEA$stat),]
gene_list_up <- deRes_GSEA$stat
names(gene_list_up) <- rownames(deRes_GSEA)
gse <- gseGO(gene_list_up,
             ont = "BP",
             OrgDb = "org.Hs.eg.db",
             keyType = "ENSEMBL",
             eps = 1e-300)
as.data.frame(gse)
#make dotplot
#make a dotplot
tiff("Plots/Control_vs_shRNA_gse_gene_ontology.tiff",width = 35,height = 21,units = 'cm',res = 300)
dotplot(gse,showCategory=30,split=".sign",
        label_format=60,font.size=7) +
  facet_grid(.~.sign)
dev.off()

tiff("Plots/Control_vs_shRNA_gse_gene_ontology_small.tiff",width = 35,height = 21,units = 'cm',res = 300)
dotplot(gse,showCategory=10,split=".sign",
        label_format=60,font.size=7) +
  facet_grid(.~.sign)
dev.off()

hs_hallmark_sets <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)
h_gene_set <- msigdbr(species = "Homo sapiens", category = "H") %>% select(gs_name,human_ensembl_gene,human_gene_symbol)

gsea_hallmark <- clusterProfiler::GSEA(gene_list_up,
                                       minGSSize = 15,
                                       maxGSSize = 500,
                                       pAdjustMethod = "BH",
                                       seed = TRUE,
                                       eps = 1e-100,
                                       nPermSimple=10000,
                                       by = "fgsea",
                                       pvalueCutoff = 1,
                                       TERM2GENE = h_gene_set)


gsea_result <- as.data.frame(gsea_hallmark)
gsea_result_summary <- gsea_result[,c(1,3,5:8)]
gsea_result_summary$adjPvalue <- ifelse(gsea_result_summary$p.adjust <= 0.1, "significant", "non-significant")
cols <- c("non-significant" ="grey","significant"="red")
#general Plot
gsea_hallmark_plot <- ggplot(gsea_result_summary, aes(reorder(ID, NES), NES, fill = adjPvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")+
  ggtitle("Control vs shRNA Hallmark Pathways")+
  theme_cowplot()
tiff("Plots/Hallmark_Gsea_bar_shRNA_vs_control.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_hallmark_plot)
dev.off()

gsea_dot <- enrichplot::dotplot(gsea_hallmark,showCategory=6,split=".sign",
                                font.size=5) +
  facet_grid(.~.sign)+
  ggtitle("Hallmark GSEA Dotplot Control vs shRNA")
tiff("Plots/Hallmark_GSEA_dot_shRNA_vs_control.tiff",width = 35,height = 21,units = 'cm',res = 300)
plot(gsea_dot)
dev.off()

## hallmark pathways gsea plots

g1 <- enrichplot::gseaplot2(gse,geneSetID = gse@result$ID[15],title = gse@result$Description[15])
tiff("Plots/Control_vs_Day5/chromatin_remodelling.tiff",width = 35,height = 21,units = 'cm',res = 300)
print(g1)
dev.off()









#plot gsea

fit3 <- gseaplot(gse,geneSetID = 34,by='all',title = gse$Description[34])
fit3
#make a pubmed trends plot
terms <- gse$Description[1:3]
pmcplot(terms,2010:2020,proportion = F)+theme_cowplot()


#make an enrichment map
x2 <- pairwise_termsim(gse)
emapplot(x2,showCategory = 20)

#make a category netplot
cnetplot(x2,categorySize="pvalue",foldChange = NULL,showCategory = 3)

# make a ridgeplot
ridgeplot(gse) + labs(x = "enrichment distribution")



# KEGG enrichment analyis
search_kegg_organism('hsa',by='kegg_code')
hs <- search_kegg_organism('Homo sapiens',by='scientific_name')
dim(hs)

# convert ensembl to interz id
down_genes_df <- as.data.frame(down_genes)
down_genes_df$entrez <- mapIds(org.Hs.eg.db,
                            keys = down_genes_df$down_genes,
                            column = "ENTREZID",
                            keytype = "ENSEMBL",
                            multiVals = "first")
genes_enterz <- mapIds(org.Hs.eg.db,
                       keys = rownames(sigOE_qValue),
                       keytype = "ENSEMBL",
                       column = "ENTREZID",
                       multiVals = "first")
genes_enterz <- unique(genes_enterz)
genes_enterz <- na.omit(genes_enterz)
kk <- enrichKEGG(gene = genes_enterz,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)
xx <- compareCluster()

head(kk)
cnetplot(kk,showCategory = 2)
dotplot(kk,showCategory=20,title="Enriched Pathways")
emapplot(pairwise_termsim(kk))

ridgeplot(kk) + labs(x = "enrichment distribution")





DEgenes=NULL
genes <- as.data.frame(as.numeric(deRes$threshold))
rownames(genes) <- rownames(deRes)
colnames(genes) <- "DEgenes"

pwf =nullp(gene.vector,"hg19","ensGene")
head(pwf)


GO.wall = goseq(pwf,"hg19","ensGene")

GO.samp=goseq(pwf,"hg19","ensGene",method = "Sampling",repcnt = 1000)
head(GO.samp)

plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.wall[,1],GO.nobias[,1]),2]),
     xlab="log10(Wallenius p-values)",ylab="log10(Sampling p-values)",
     xlim=c(-3,0))
abline(0,1,col=3,lty=2)



GO.nobias=goseq(pwf,"hg19","ensGene",method="Hypergeometric")
head(GO.nobias)



enriched.GO=GO.nobias$category[GO.nobias$over_represented_pvalue<0.01]
enriched.GO=GO.nobias$category[p.adjust(GO.nobias$over_represented_pvalue,
                                       method="BH")<0.05]
head(enriched.GO)
length(enriched.GO)


library(GO.db)

capture.output(for(go in enriched.GO[1:32]) { print(GOTERM[[go]])
  cat("--------------------------------------\n")
}
, file="SigGo.txt")

# Get the mapping from ENSEMBL 2 Entrez
en2eg=as.list(org.Hs.egENSEMBL2EG)
# Get the mapping from Entrez 2 KEGG
eg2kegg=as.list(org.Hs.egPATH)
# Define a function which gets all unique KEGG IDs
# associated with a set of Entrez IDs
grepKEGG=function(id,mapkeys){unique(unlist(mapkeys[id],use.names=FALSE))}
# Apply this function to every entry in the mapping from
# ENSEMBL 2 Entrez to combine the two maps
kegg=lapply(en2eg,grepKEGG,eg2kegg)
head(kegg)

KEGG = goseq(pwf,gene2cat = kegg)

####

ego <- clusterProfiler::enrichGO(gene = rownames(sigOE_qValue),
                universe = rownames(deRes),
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.05,
                readable = TRUE)

cluster_summary <- data.frame(ego)

dotplot(ego,showCategory=50)
x2 <- enrichplot::pairwise_termsim(ego) 
emapplot(x2)
emapplot(ego,showCategory = 50)

ego <- enrichKEGG(ego)

cnetplot(ego,
         categorySize="pvalue",
         showCategory = 50,
         foldChange= sigOE_qValue$log2FoldChange,
         vertex.label.font=6)

## Enrichment analysis for all control vs shRNA,


