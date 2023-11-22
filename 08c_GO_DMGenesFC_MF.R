#This code performs the GO term analysis of the DM genes
#Uses the list generated in 06_featureCountsGeneLists.R
#Author: Lea Faivre
#Date: 231017

# Libraries ---------------------------------------------------------------

library(dplyr)
library(readxl)
library(tidyverse)
library(ggplot2)
library(rio)
library(ggvenn)
library(cowplot)
library(xlsx)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(topGO)
library(org.At.tair.db)

# Data import -------------------------------------------------------------
DMGenes <- import_list("Data/DMGenes_FC.xlsx", setclass = "tbl")

Araport <- TxDb.Athaliana.BioMart.plantsmart28
Genes <- GenomicFeatures::genes(Araport,columns = c("GENEID","TXTYPE"))
Genesdf <- as.data.frame(Genes)

# GO analysis biological function -------------------------------------------------------------

#Making a vector of interesting genes: all genes, with 0 for non DM and 1 for DM
GeneList <- lapply(DMGenes, function(x){
  a <- factor(as.integer(Genesdf$GENEID %in% x$Gene))
  names(a) <- Genesdf$GENEID
  a
})

GOTest <- lapply(GeneList, function(a){
  new("topGOdata", ontology = "MF", allGenes = a, 
      annot = annFUN.org, mapping = "org.At.tair.db", geneSelectionFun = function(x)(x))
})


GOTable <- lapply(GOTest, function(x){
  Res <- runTest(x, algorithm = "weight01", statistic = "fisher")
  GenTable(x, raw.p.value = Res,
           topNodes = length(Res@score), numChar = 60)
})

FormattedTable <- lapply(GOTable, function(x){
  x %>% 
    mutate(raw.p.value = as.numeric(raw.p.value),
           logPvalue = log10(1/raw.p.value),
           FoldEnrichment = Significant/Expected) %>%
    filter(raw.p.value <= 0.01)
})


for (i in names(FormattedTable)) {
  write.xlsx(x = as.data.frame(FormattedTable[[i]]),
  file = "Results/GO_DMGenes_FC_MF.xlsx",
  sheetName = paste(i, sep = ""),
  row.names = FALSE,
  append = T)
}