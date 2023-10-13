#This code plots the Venn diagrams of the DM genes
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