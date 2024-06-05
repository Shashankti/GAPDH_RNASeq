library(plyr)
library(LSD)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(stringr)
library(topGO)
library(genefilter)
library(biomaRt)
library(dplyr)
library(EDASeq)
library(fdrtool)
library(org.Hs.eg.db)

overallBaseMean <- as.matrix(deRes[,"baseMean",drop=F])
anSig <- 
backG <- genefinder(overallBaseMean,rownames(sigOE_down),10,method = "manhattan")
## get identified similar genes
backG <- rownames(overallBaseMean)[as.vector(sapply(backG, function(x)x$indices))]
## remove DE genes from background
backG <- setdiff(backG, rownames(sigOE_down))
## number of genes in the background
length(backG)

multidensity(list(
  all = log2(deRes$baseMean),
  fore = log2(deRes[rownames(sigOE_down),"baseMean"]),
  back = log2(deRes[backG,"baseMean"])),
  xlab ="log2 mean counts", main = "Matching for enrichment analysis"
)
ggplot(deRes,aes(x=log2(baseMean)))+geom_density()

onts = c("MF","BP","CC")
geneIDs <- rownames(overallBaseMean)
inUni <- geneIDs %in% c(rownames(sigOE_down),backG)
inSel <- geneIDs %in% rownames(sigOE_down)

alg <- factor(as.integer(inSel[inUni]))
names(alg) <- geneIDs[inUni]

tab = as.list(onts)
names(tab) = onts

for(i in 1:3){
  ## prepare data
  tgd <- new( "topGOdata", ontology=onts[i], allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )
  ## run tests
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  ## look at results
  tab[[i]] <- GenTable( tgd, Fisher.elim = resultTopGO.elim,
                        Fisher.classic = resultTopGO.classic,
                        orderBy = "Fisher.classic" , topNodes = 200)
}

topGOResults <- rbind.fill(tab)

#######################################333

#using GSEA

