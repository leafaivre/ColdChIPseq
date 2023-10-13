#This code performs the GO term analysis of the DM genes
#Uses the list generated in 03_GenerationGeneLists
#Author: Lea Faivre
#Date: 231013

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
DMGenes <- import_list("Data/DMGenes_Full.xlsx", setclass = "tbl")

UniqueDM <- lapply(DMGenes, function(x){
  y <- x %>%
    dplyr::select(Gene)
  y <- unique(y)
})

Araport <- TxDb.Athaliana.BioMart.plantsmart28
Genes <- GenomicFeatures::genes(Araport,columns = c("GENEID","TXTYPE"))
Genesdf <- as.data.frame(Genes)

# GO analysis -------------------------------------------------------------

DiffGenes <- list(K4hU = anti_join(UniqueDM$K4_3h_Gain, UniqueDM$K4_3h_Loss, by = "Gene"),
                  K4hD = anti_join(UniqueDM$K4_3h_Loss, UniqueDM$K4_3h_Gain, by = "Gene"),
                  K4dU = anti_join(UniqueDM$K4_3d_Gain, UniqueDM$K4_3d_Loss, by = "Gene"),
                  K4dD = anti_join(UniqueDM$K4_3d_Loss, UniqueDM$K4_3d_Gain, by = "Gene"),
                  K27hU = anti_join(UniqueDM$K27_3h_Gain, UniqueDM$K27_3h_Loss, by = "Gene"),
                  K27hD = anti_join(UniqueDM$K27_3h_Loss, UniqueDM$K27_3h_Gain, by = "Gene"),
                  K27dU = anti_join(UniqueDM$K27_3d_Gain, UniqueDM$K27_3d_Loss, by = "Gene"),
                  K27dD = anti_join(UniqueDM$K27_3d_Loss, UniqueDM$K27_3d_Gain, by = "Gene"))

#Making a vector of interesting genes: all genes, with 0 for non DM and 1 for DM
GeneList <- lapply(DiffGenes, function(x){
  a <- factor(as.integer(Genesdf$GENEID %in% x$Gene))
  names(a) <- Genesdf$GENEID
  a
})

GOTest <- lapply(GeneList, function(a){
  new("topGOdata", ontology = "BP", allGenes = a, 
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


for (i in names(SigGO)) {
  write.xlsx(x = as.data.frame(SigGO[[i]]),
             file = "Results/GO_DMGenes_woUPDOWN_0.01.xlsx",
             sheetName = paste(i, sep = ""),
             row.names = FALSE,
             append = T)
}


# Extracting Genes related to a Term --------------------------------------

GO <- lapply(GOTest, genesInTerm)

TAIR10genes <- read_excel("InputFiles/TAIR10genes.xlsx", 
                          sheet = "Name", col_types = c("text", "text"))

ColdRelated <- lapply(DiffGenes, function(x) {
  x %>%
    filter(Gene %in% unlist(allGO$K4hU["GO:0009409"])) %>%
    left_join(TAIR10genes, by = "Gene")
  
})



