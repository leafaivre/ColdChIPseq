#This code performs the GO term analysis of the DE genes
#It uses the DE gene lists generated in 10_RNAseqAnalysis
#Author: Lea Faivre
#Date: 231031

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
DEGenes <- import_list("Data/DEGenes.xlsx", setclass = "tbl")

DE <- list("h_UP" = DEGenes$Nvs3h %>%
             filter(log2FoldChange > 1),
           "h_DOWN" = DEGenes$Nvs3h %>%
             filter(log2FoldChange < -1),
           "d_UP" = DEGenes$Nvs3d %>%
             filter(log2FoldChange > 1),
           "d_DOWN" = DEGenes$Nvs3d %>%
             filter(log2FoldChange < -1))

Araport <- TxDb.Athaliana.BioMart.plantsmart28
Genes <- GenomicFeatures::genes(Araport,columns = c("GENEID","TXTYPE"))
Genesdf <- as.data.frame(Genes)

# GO analysis biological function -------------------------------------------------------------

#Making a vector of interesting genes: all genes, with 0 for non DM and 1 for DM
GeneList <- lapply(DE, function(x){
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
  file = "Results/GO_DEGenes_FC_0.01.xlsx",
  sheetName = paste(i, sep = ""),
  row.names = FALSE,
  append = T)
}

# Graphs ------------------------------------------------------------------

Top20 <- lapply(FormattedTable, function(x){
  x %>%
    filter(FoldEnrichment > 1) %>%
    slice_max(logPvalue, n = 20, with_ties = FALSE) %>%
    mutate(GOTerm = paste(Term, "(", GO.ID, ")"))
})

ggplot(Top20$h_UP, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                             fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("skyblue", "plum4", "skyblue","skyblue", "skyblue",
                               "skyblue", "skyblue", "skyblue","skyblue", "dimgray",
                               "dimgray", "royalblue4", "skyblue","dimgray", "skyblue",
                               "skyblue", "dimgray", "dimgray","dimgray", "skyblue")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes upregulated after 3h of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO/3h_DE_Up.tiff", units = "in", width = 8, height = 4, dpi = 300, compression = 'lzw')

ggplot(Top20$h_DOWN, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                             fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("skyblue", "dimgray", "dimgray","dimgray", "skyblue",
                               "darkseagreen", "skyblue", "dimgray","dimgray", "dimgray",
                               "dimgray", "dimgray", "dimgray","dimgray", "dimgray",
                               "skyblue", "darkseagreen", "plum4","skyblue", "dimgray")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes downregulated after 3h of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO/3h_DE_Down.tiff", units = "in", width = 8, height = 3, dpi = 300, compression = 'lzw')

ggplot(Top20$d_UP, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                             fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("skyblue", "royalblue4", "skyblue", "skyblue", "skyblue",
                               "skyblue", "skyblue", "skyblue", "skyblue", "dimgray",
                               "skyblue", "skyblue", "dimgray", "dimgray", "skyblue",
                               "dimgray", "dimgray", "skyblue", "dimgray", "skyblue")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes upregulated after 3d of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO/3d_DE_Up.tiff", units = "in", width = 8, height = 4, dpi = 300, compression = 'lzw')

ggplot(Top20$d_DOWN, aes(x = reorder(GOTerm, logPvalue), y = logPvalue, 
                             fill = reorder(GOTerm, -logPvalue))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("skyblue", "dimgray", "dimgray", "coral3", "dimgray",
                               "skyblue", "dimgray", "skyblue", "dimgray", "skyblue",
                               "skyblue", "dimgray", "skyblue", "skyblue", "dimgray",
                               "skyblue", "coral3", "darkseagreen", "skyblue", "darkseagreen")) +
  labs(x = "", y = "log10(1/P-value)", title = "Genes downregulated after 3d of cold exposure") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/GO/3d_DE_Down.tiff", units = "in", width = 8, height = 4, dpi = 300, compression = 'lzw')
