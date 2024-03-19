#This code plots the UpsetR diagrams of the DM genes
#Uses the list generated in 03_GenerationGeneLists
#Author: Lea Faivre
#Date: 240319

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
library(UpSetR)

# Data import -------------------------------------------------------------
DMGenes <- import_list("Data/DMGenes_FC.xlsx", setclass = "tbl")

UniqueDM <- lapply(DMGenes, function(x){
  y <- x %>%
    dplyr::select(Gene)
  y <- unique(y)
})

# Overlaps ----------------------------------------------------------------
Ph <- list("H3K4me3 Gain" = UniqueDM$K4_3h_Gain$Gene,  
           "H3K4me3 Loss" = UniqueDM$K4_3h_Loss$Gene, 
           "H3K27me3 Gain" = UniqueDM$K27_3h_Gain$Gene, 
           "H3K27me3 Loss" = UniqueDM$K27_3h_Loss$Gene)

tiff(file = "Plots/Upset_K4K27_3h.tiff", width = 5.5, height = 5, units = "in", res = 300)

upset(fromList(Ph), order.by = "freq", 
      sets.bar.color = c("darkgreen","#749140", "darkred", "#cd6155"),
      mainbar.y.label = "Differentially methylated genes", sets.x.label = "DM Genes",
      text.scale = c(1.8, 1.8, 1.5, 1.3, 1.8, 2), point.size = 3, line.size = 1,  mb.ratio = c(0.7, 0.3),
      main.bar.color = "black", matrix.color = "black",
      queries = list(list(query = intersects, 
                          params = list("H3K4me3 Gain", "H3K27me3 Gain"), 
                          color = "darkgoldenrod", active = T), 
                     list(query = intersects, 
                          params = list("H3K4me3 Loss", "H3K27me3 Loss"), 
                          color = "darkgoldenrod", active = T), 
                     list(query = intersects, 
                          params = list("H3K4me3 Gain", "H3K27me3 Loss"), 
                          color = "#5a427e", active = T),
                     list(query = intersects, 
                          params = list("H3K4me3 Loss", "H3K27me3 Gain"), 
                          color = "#5a427e", active = T)))
dev.off()


Pd <- list("H3K4me3 Gain" = UniqueDM$K4_3d_Gain$Gene,  
           "H3K4me3 Loss" = UniqueDM$K4_3d_Loss$Gene, 
           "H3K27me3 Gain" = UniqueDM$K27_3d_Gain$Gene, 
           "H3K27me3 Loss" = UniqueDM$K27_3d_Loss$Gene)

tiff(file = "Plots/Upset_K4K27_3d.tiff", width = 5.5, height = 5, units = "in", res = 300)

upset(fromList(Pd), order.by = "freq", 
      sets = c("H3K4me3 Gain", "H3K4me3 Loss", "H3K27me3 Gain", "H3K27me3 Loss"),
      sets.bar.color = c("darkgreen","#749140", "darkred", "#cd6155"),
      keep.order = T, 
      mainbar.y.label = "Differentially methylated genes", sets.x.label = "DM Genes",
      text.scale = c(1.8, 1.8, 1.5, 1.3, 1.8, 2), point.size = 3, line.size = 1,  mb.ratio = c(0.7, 0.3),
      main.bar.color = "black", matrix.color = "black",
      queries = list(list(query = intersects, 
                          params = list("H3K4me3 Gain", "H3K27me3 Gain"), 
                          color = "darkgoldenrod", active = T), 
                     list(query = intersects, 
                          params = list("H3K4me3 Loss", "H3K27me3 Loss"), 
                          color = "darkgoldenrod", active = T), 
                     list(query = intersects, 
                          params = list("H3K4me3 Gain", "H3K27me3 Loss"), 
                          color = "#5a427e", active = T),
                     list(query = intersects, 
                          params = list("H3K4me3 Loss", "H3K27me3 Gain"), 
                          color = "#5a427e", active = T)))
dev.off()
