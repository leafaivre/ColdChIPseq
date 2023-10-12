#This code plots the Venn diagrams and GO terms of the DM genes
#Uses the list generated in 03_GenerationGeneLists
#Author: Lea Faivre
#Date: 231012

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
Genes <- genes(Araport,columns = c("GENEID","TXTYPE"))
Genesdf <- as.data.frame(Genes)

# Persistence -------------------------------------------------------------
K4_gain <- list("3h" = UniqueDM$K4_3h_Gain$Gene, "3d" = UniqueDM$K4_3d_Gain$Gene)
K4_loss <- list("3h" = UniqueDM$K4_3h_Loss$Gene, "3d" = UniqueDM$K4_3d_Loss$Gene)
K27_gain <- list("3h" = UniqueDM$K27_3h_Gain$Gene, "3d" = UniqueDM$K27_3d_Gain$Gene)
K27_loss <- list("3h" = UniqueDM$K27_3h_Loss$Gene, "3d" = UniqueDM$K27_3d_Loss$Gene)

K4_L <- ggvenn(K4_loss, fill_color = c("#cd6155","#5f0000"), 
               stroke_linetype = "blank", 
               show_percentage = FALSE, fill_alpha = 0.7)

K4_G <- ggvenn(K4_gain, fill_color = c("#749140","#196F3D"), 
               stroke_linetype = "blank", 
               show_percentage = FALSE, fill_alpha = 0.7)

K27_L <- ggvenn(K27_loss, fill_color = c("#cd6155","#5f0000"), 
                stroke_linetype = "blank", 
                show_percentage = FALSE, fill_alpha = 0.7)

K27_G <- ggvenn(K27_gain, fill_color = c("#749140","#196F3D"), 
                stroke_linetype = "blank",
                show_percentage = FALSE, fill_alpha = 0.7)

grid.newpage()
plot_grid(K4_L, K4_G, K27_L, K27_G, 
          labels = c("K4me3 Loss", "K4me3 gain", "K27me3 Loss", "K27me3 Gain"), 
          label_x = 0.5, hjust = 0.5)
#ggsave("Plots/PersistenceDMG.tiff", units = "in", width = 7, height = 7, 
       #dpi = 300, compression = 'lzw')


# Overlaps ----------------------------------------------------------------
Ph <- list("H3K4me3 Gain" = UniqueDM$K4_3h_Gain$Gene,  
           "H3K4me3 Loss" = UniqueDM$K4_3h_Loss$Gene, 
           "H3K27me3 Gain" = UniqueDM$K27_3h_Gain$Gene, 
           "H3K27me3 Loss" = UniqueDM$K27_3h_Loss$Gene)

ggvenn(Ph, fill_color = c("#749140","sandybrown", "#cd6155", "aquamarine3"), 
       stroke_linetype = "blank", show_percentage = FALSE, 
       fill_alpha = 0.7, text_size = 4, set_name_size = 4)
#ggsave("Plots/K4K27_3h.tiff", units = "in", width = 10, height = 5, dpi = 300, 
       #compression = 'lzw')


Pd <- list("H3K4me3 Gain" = UniqueDM$K4_3d_Gain$Gene,  
           "H3K4me3 Loss" = UniqueDM$K4_3d_Loss$Gene, 
           "H3K27me3 Gain" = UniqueDM$K27_3d_Gain$Gene, 
           "H3K27me3 Loss" = UniqueDM$K27_3d_Loss$Gene)

ggvenn(Pd, fill_color = c("#749140","sandybrown", "#cd6155", "aquamarine3"), 
       stroke_linetype = "blank", show_percentage = FALSE, 
       fill_alpha = 0.7, text_size = 4, set_name_size = 4)
#ggsave("Plots/K4K27_3d.tiff", units = "in", width = 10, height = 5, dpi = 300, 
       #compression = 'lzw')


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

