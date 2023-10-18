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

# GO analysis -------------------------------------------------------------

#Making a vector of interesting genes: all genes, with 0 for non DM and 1 for DM
GeneList <- lapply(DMGenes, function(x){
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


for (i in names(FormattedTable)) {
  write.xlsx(x = as.data.frame(FormattedTable[[i]]),
             file = "Results/GO_DMGenes_FC_0.01.xlsx",
             sheetName = paste(i, sep = ""),
             row.names = FALSE,
             append = T)
}


# Extracting Genes related to a Term --------------------------------------

allGO <- genesInTerm(GOTest$K4_3h_Gain)


  lapply(GOTest, genesInTerm)

Interesting <- lapply(DMGenes, function(x) {
  x %>%
    filter(Gene %in% unlist(allGO["GO:0009266"])) 
})


# Graphs ------------------------------------------------------------------

Top20 <- lapply(FormattedTable, function(x){
  x %>%
    filter(Significant >= 2) %>%
    filter(FoldEnrichment > 1) %>%
    slice_max(logPvalue, n = 20, with_ties = FALSE) %>%
    mutate(GOTerm = paste(Term, "(", GO.ID, ")"))
  
  
})

ggplot(Top20$K4_3h_Gain, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                           fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("skyblue", "dimgray", "dimgray","skyblue", "royalblue4",
                               "dimgray", "skyblue", "dimgray","coral3", "skyblue",
                               "dimgray", "royalblue4", "dimgray","skyblue", "dimgray",
                               "skyblue", "dimgray", "skyblue","dimgray", "plum4")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes gaining H3K4me3 after 3h of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO/3h_K4_Gain.tiff", units = "in", width = 8, height = 4, dpi = 300, compression = 'lzw')

ggplot(Top20$K4_3h_Loss, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                             fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("dimgray", "darkseagreen", "plum4", "dimgray", "dimgray",
                               "darkseagreen", "coral3", "skyblue", "plum4", "dimgray",
                               "royalblue4", "skyblue", "dimgray", "dimgray", "darkseagreen",
                               "dimgray", "dimgray", "plum4", "dimgray", "dimgray")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes losing H3K4me3 after 3h of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO/3h_K4_Loss.tiff", units = "in", width = 8, height = 4, dpi = 300, compression = 'lzw')

ggplot(Top20$K4_3d_Gain, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                             fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("dimgray", "dimgray", "royalblue4", "skyblue", "plum4",
                               "skyblue", "skyblue", "skyblue", "skyblue", "dimgray",
                               "royalblue4", "skyblue", "dimgray", "royalblue4", "skyblue",
                               "coral3", "plum4", "skyblue", "dimgray", "coral3")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes gaining H3K4me3 after 3d of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO/3d_K4_Gain.tiff", units = "in", width = 8, height = 4, dpi = 300, compression = 'lzw')

ggplot(Top20$K4_3d_Loss, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                             fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("coral3", "darkseagreen", "coral3", "dimgray", "skyblue",
                               "royalblue4", "goldenrod3", "skyblue", "coral3", "dimgray",
                               "darkseagreen", "skyblue", "dimgray", "coral3", "skyblue",
                               "goldenrod3", "dimgray", "darkseagreen", "dimgray", "skyblue")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes losing H3K4me3 after 3d of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO/3d_K4_Loss.tiff", units = "in", width = 8, height = 4, dpi = 300, compression = 'lzw')

ggplot(Top20$K27_3h_Gain, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                             fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("plum4", "dimgray", "dimgray", "plum4", "plum4",
                               "dimgray", "dimgray", "darkseagreen", "darkseagreen", "darkseagreen",
                               "darkseagreen", "darkseagreen", "darkseagreen", "darkseagreen", "darkseagreen",
                               "goldenrod3", "dimgray", "darkseagreen", "dimgray", "dimgray")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes gaining H3K27me3 after 3h of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO/3h_K27_Gain.tiff", units = "in", width = 8, height = 4, dpi = 300, compression = 'lzw')

ggplot(Top20$K27_3h_Loss, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                              fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("dimgray", "dimgray", "dimgray", "dimgray", "goldenrod3",
                               "dimgray", "dimgray", "plum4", "skyblue", "dimgray", "dimgray")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes losing H3K27me3 after 3h of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO/3h_K27_Loss.tiff", units = "in", width = 8, height = 2.8, dpi = 300, compression = 'lzw')

ggplot(Top20$K27_3d_Gain, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                              fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("plum4", "plum4", "dimgray", "plum4", "dimgray",
                               "darkseagreen", "dimgray", "dimgray", "darkseagreen", "dimgray",
                               "darkseagreen", "dimgray", "dimgray", "goldenrod3", "darkseagreen",
                               "darkseagreen", "dimgray", "darkseagreen", "skyblue", "darkseagreen")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes gaining H3K27me3 after 3d of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO/3d_K27_Gain.tiff", units = "in", width = 8, height = 4, dpi = 300, compression = 'lzw')

ggplot(Top20$K27_3d_Loss, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                              fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("dimgray", "dimgray", "dimgray", "dimgray", "dimgray",
                               "dimgray", "dimgray", "skyblue", "dimgray", "dimgray",
                               "dimgray", "dimgray", "skyblue", "skyblue", "dimgray")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes losing H3K27me3 after 3d of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
ggsave("Plots/GO/3d_K27_Loss.tiff", units = "in", width = 8, height = 3, dpi = 300, compression = 'lzw')
