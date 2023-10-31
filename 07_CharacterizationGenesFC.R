#This code plots the Venn diagrams of the DM genes
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

Pie <- import_list("Data/Persistence_Pie.xlsx", setclass = "tbl")

# Persistence -------------------------------------------------------------
K4_gain <- list("3h" = DMGenes$K4_3h_Gain$Gene, "3d" = DMGenes$K4_3d_Gain$Gene)
K4_loss <- list("3h" = DMGenes$K4_3h_Loss$Gene, "3d" = DMGenes$K4_3d_Loss$Gene)
K27_gain <- list("3h" = DMGenes$K27_3h_Gain$Gene, "3d" = DMGenes$K27_3d_Gain$Gene)
K27_loss <- list("3h" = DMGenes$K27_3h_Loss$Gene, "3d" = DMGenes$K27_3d_Loss$Gene)

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


# Persistence Pie Charts ----------------------------------------------------------------------

Pie <- plyr::ldply(Pie, data.frame, .id = "List") %>%
  group_by(List) %>%
  mutate(Percentage = Nb / sum(Nb),
         Cat = as.factor(Cat),
         Cat = fct_relevel(Cat, "3h only", "3d only", "3h & 3d"),
         Label = paste(round(Percentage*100, digits = 1), "%", sep = ""))
  
ggplot(Pie, aes(x = "", y = Percentage, fill = Cat)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("#92C5DE","#2166AC", "darkgoldenrod")) +
  geom_label(aes(label = Label, x = 1.15),
             position = position_stack(vjust = 0.5),
             label.size = NA) +
  coord_polar("y", start = 0) +
  facet_wrap(~List, ncol = 2) +
  theme_void() +
  theme(legend.position = "none")
ggsave("Plots/PersistencePie.tiff", units = "in", width = 6, height = 6, 
       dpi = 300, compression = 'lzw')

# Overlaps ----------------------------------------------------------------
Ph <- list("H3K4me3 Gain" = DMGenes$K4_3h_Gain$Gene,  
           "H3K4me3 Loss" = DMGenes$K4_3h_Loss$Gene, 
           "H3K27me3 Gain" = DMGenes$K27_3h_Gain$Gene, 
           "H3K27me3 Loss" = DMGenes$K27_3h_Loss$Gene)

ggvenn(Ph, fill_color = c("#749140","sandybrown", "#cd6155", "aquamarine3"), 
       stroke_linetype = "blank", show_percentage = FALSE, 
       fill_alpha = 0.7, text_size = 4, set_name_size = 4)
ggsave("Plots/K4K27_3h.tiff", units = "in", width = 3, height = 5, dpi = 300, compression = 'lzw')


Pd <- list("H3K4me3 Gain" = DMGenes$K4_3d_Gain$Gene,  
           "H3K4me3 Loss" = DMGenes$K4_3d_Loss$Gene, 
           "H3K27me3 Gain" = DMGenes$K27_3d_Gain$Gene, 
           "H3K27me3 Loss" = DMGenes$K27_3d_Loss$Gene)

ggvenn(Pd, fill_color = c("#465726","#da6d10", "#893228", "#37a581"), 
       stroke_linetype = "blank", show_percentage = FALSE, 
       fill_alpha = 0.7, text_size = 4, set_name_size = 4)
ggsave("Plots/K4K27_3d.tiff", units = "in", width = 3, height = 5, dpi = 300, compression = 'lzw')
