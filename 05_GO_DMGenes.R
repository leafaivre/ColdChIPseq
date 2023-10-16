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
DMGenes <- import_list("Data/Test/DMGenes_Full.xlsx", setclass = "tbl")

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

SigGO <- lapply(GOTable, function(x){
  x %>% 
    mutate(raw.p.value = as.numeric(raw.p.value),
           logPvalue = log10(1/raw.p.value),
           FoldEnrichment = Significant/Expected) %>%
    filter(raw.p.value <= 0.01)
})


for (i in names(SigGO)) {
  write.xlsx(x = as.data.frame(SigGO[[i]]), 
             file = "Results/New_GO_DMGenes.xlsx", sheetName = paste(i, sep = ""),
             row.names = FALSE, append = T)}


# Extracting Genes related to a Term --------------------------------------

allGO <- lapply(GOTest, genesInTerm)

TAIR10genes <- read_excel("InputFiles/TAIR10genes.xlsx", 
                          sheet = "Name", col_types = c("text", "text"))

ColdRelated <- lapply(DiffGenes, function(x) {
  x %>%
    filter(Gene %in% unlist(allGO$K4hU["GO:0009409"])) %>%
    left_join(TAIR10genes, by = "Gene")
  
})


# Plotting ----------------------------------------------------------------
Top20 <- lapply(SigGO, function(x){
  x %>%
    filter(Significant >= 2) %>%
    filter(FoldEnrichment >= 1.5) %>%
    slice_max(logPvalue, n = 20, with_ties = FALSE) %>%
    mutate(GOTerm = paste(Term, "(", GO.ID, ")"))
})

ggplot(Top20$K4hU, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                           fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("#92C5DE", "#92C5DE", "#2166AC", "dimgray", "#92C5DE", 
                               "dimgray", "coral3", "#2166AC","#92C5DE", "dimgray",
                               "dimgray", "coral3",  "dimgray", "coral3", "#92C5DE",
                              "dimgray", "dimgray","dimgray", "#92C5DE", "#92C5DE")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes gaining H3K4me3 after 3h") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO_3h_K4_Gain.tiff", units = "in", width = 8, height = 4, dpi = 300, compression = 'lzw')

ggplot(Top20$K4hD, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                       fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("dimgray", "darkseagreen", "dimgray", "coral3", "darkseagreen", 
                               "goldenrod2", "#92C5DE", "#2166AC","goldenrod2", "dimgray",
                               "dimgray", "coral3",  "dimgray", "dimgray", "#92C5DE",
                               "goldenrod3", "coral3","goldenrod3", "goldenrod3", "dimgray")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes losing H3K4me3 after 3h") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO_3h_K4_Loss.tiff", units = "in", width = 8, height = 4, dpi = 300, compression = 'lzw')
