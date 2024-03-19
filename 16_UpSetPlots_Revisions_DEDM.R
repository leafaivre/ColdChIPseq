#This code plots the Upset diagram between DM and DE genes
#It uses the DE gene lists generated in 10_RNAseqAnalysis and 06_featureCountsGeneLists
#Author: Lea Faivre
#Date: 240319


# Libraries -----------------------------------------------------------------------------------

library(readxl)
library(tidyverse)
library(RColorBrewer)
library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rio)
library(ggvenn)
library(UpSetR)

# Data import ---------------------------------------------------------------------------------

DM_FC <- import_list("Data/FoldChanges_allGenes.xlsx", setclass = "tbl")

DM_FC <- lapply(DM_FC, function(x){
  x %>% 
    select(Gene, log2FCNvs3h, log2FCNvs3d)
})

DM_FC_form <- list(
  K4 = DM_FC$K4 %>%
    dplyr::rename(K4_FC_3h = log2FCNvs3h, K4_FC_3d = log2FCNvs3d),
  K27 = DM_FC$K27 %>%
    dplyr::rename(K27_FC_3h = log2FCNvs3h, K27_FC_3d = log2FCNvs3d))

DE_FC <- list(
  h = read_excel("Data/RNAseq_Nvs3h.xlsx") %>%
    select(Gene, log2FoldChange, padj) %>%
    dplyr::rename(DE_FC_3h = log2FoldChange, padj_3h = padj),
  d = read_excel("Data/RNAseq_Nvs3d.xlsx") %>%
    select(Gene, log2FoldChange, padj) %>%
    dplyr::rename(DE_FC_3d = log2FoldChange, padj_3d = padj))

DMGenes <- import_list("Data/DMGenes_FC.xlsx", setclass = "tbl")

DEGenes <- import_list("Data/DEGenes.xlsx", setclass = "tbl")

DE <- list("h_UP" = DEGenes$Nvs3h %>%
             filter(log2FoldChange > 1),
           "h_DOWN" = DEGenes$Nvs3h %>%
             filter(log2FoldChange < -1),
           "d_UP" = DEGenes$Nvs3d %>%
             filter(log2FoldChange > 1),
           "d_DOWN" = DEGenes$Nvs3d %>%
             filter(log2FoldChange < -1))

# 3h UP -------------------------------------------------------------------------------
h_UP <- list("Gaining H3K4me3" = DMGenes$K4_3h_Gain$Gene, 
             "Losing H3K27me3" = DMGenes$K27_3h_Loss$Gene,
             "Upregulated" = DE$h_UP$Gene)

tiff(file = "Plots/UpsetDEDM/3h_UP.tiff", width = 5.5, height = 5, units = "in", res = 300)

upset(fromList(h_UP), order.by = "freq", 
      sets.bar.color = c("#aa73aa", "darkgreen", "#cd6155"),
      sets = c("Upregulated", "Gaining H3K4me3", "Losing H3K27me3"),
      keep.order = T,
      mainbar.y.label = "DE and DM genes", sets.x.label = "DM/DE Genes",
      text.scale = c(1.8, 1.8, 1.5, 1.3, 1.8, 2), point.size = 3, line.size = 1,  mb.ratio = c(0.7, 0.3),
      main.bar.color = "black", matrix.color = "black",
      queries = list(list(query = intersects, 
                          params = list("Gaining H3K4me3", "Upregulated"), 
                          color = "darkgreen", active = T), 
                     list(query = intersects, 
                          params = list("Losing H3K27me3", "Upregulated"), 
                          color = "#cd6155", active = T), 
                     list(query = intersects, 
                          params = list("Gaining H3K4me3", "Losing H3K27me3", "Upregulated"), 
                          color = "#aa73aa", active = T)))
dev.off()


# 3h DOWN -------------------------------------------------------------------------------------
h_DOWN <- list("Losing H3K4me3" = DMGenes$K4_3h_Loss$Gene, 
               "Gaining H3K27me3" = DMGenes$K27_3h_Gain$Gene,
               "Downregulated" = DE$h_DOWN$Gene)

tiff(file = "Plots/UpsetDEDM/3h_DOWN.tiff", width = 5.5, height = 5, units = "in", res = 300)

upset(fromList(h_DOWN), order.by = "freq",
      sets = c("Downregulated", "Losing H3K4me3", "Gaining H3K27me3"),
      keep.order = T,
      sets.bar.color = c("powderblue", "#749140", "darkred"),
      mainbar.y.label = "DE and DM genes", sets.x.label = "DM/DE Genes",
      text.scale = c(1.8, 1.8, 1.5, 1.3, 1.8, 2), point.size = 3, line.size = 1,  mb.ratio = c(0.7, 0.3),
      main.bar.color = "black", matrix.color = "black",
      queries = list(list(query = intersects, 
                          params = list("Losing H3K4me3", "Downregulated"), 
                          color = "#749140", active = T), 
                     list(query = intersects, 
                          params = list("Gaining H3K27me3", "Downregulated"), 
                          color = "darkred", active = T)))

dev.off()


# 3d UP ---------------------------------------------------------------------------------------
d_UP <- list("Gaining H3K4me3" = DMGenes$K4_3d_Gain$Gene, 
             "Losing H3K27me3" = DMGenes$K27_3d_Loss$Gene,
             "Upregulated" = DE$d_UP$Gene)

tiff(file = "Plots/UpsetDEDM/3d_UP.tiff", width = 5.5, height = 5, units = "in", res = 300)

upset(fromList(d_UP), order.by = "freq", 
      sets.bar.color = c("#653e65", "darkgreen", "#cd6155"),
      sets = c("Upregulated", "Gaining H3K4me3", "Losing H3K27me3"),
      keep.order = T,
      mainbar.y.label = "DE and DM genes", sets.x.label = "DM/DE Genes",
      text.scale = c(1.8, 1.8, 1.5, 1.3, 1.8, 2), point.size = 3, line.size = 1,  mb.ratio = c(0.7, 0.3),
      main.bar.color = "black", matrix.color = "black",
      queries = list(list(query = intersects, 
                          params = list("Gaining H3K4me3", "Upregulated"), 
                          color = "darkgreen", active = T), 
                     list(query = intersects, 
                          params = list("Losing H3K27me3", "Upregulated"), 
                          color = "#cd6155", active = T), 
                     list(query = intersects, 
                          params = list("Gaining H3K4me3", "Losing H3K27me3", "Upregulated"), 
                          color = "#653e65", active = T)))
dev.off()


# 3d DOWN -------------------------------------------------------------------------------------
d_DOWN <- list("Losing H3K4me3" = DMGenes$K4_3d_Loss$Gene, 
               "Gaining H3K27me3" = DMGenes$K27_3d_Gain$Gene,
               "Downregulated" = DE$d_DOWN$Gene)

tiff(file = "Plots/UpsetDEDM/3d_DOWN.tiff", width = 5.5, height = 5, units = "in", res = 300)

upset(fromList(d_DOWN), order.by = "freq",
      sets = c("Downregulated", "Losing H3K4me3", "Gaining H3K27me3"),
      keep.order = T,
      sets.bar.color = c("#27727b", "#749140", "darkred"),
      mainbar.y.label = "DE and DM genes", sets.x.label = "DM/DE Genes",
      text.scale = c(1.8, 1.8, 1.5, 1.3, 1.8, 2), point.size = 3, line.size = 1,  mb.ratio = c(0.7, 0.3),
      main.bar.color = "black", matrix.color = "black",
      queries = list(list(query = intersects, 
                          params = list("Losing H3K4me3", "Downregulated"), 
                          color = "#749140", active = T), 
                     list(query = intersects, 
                          params = list("Gaining H3K27me3", "Downregulated"), 
                          color = "darkred", active = T)))

dev.off()
