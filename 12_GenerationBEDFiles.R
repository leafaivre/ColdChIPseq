#This code generates the BED files needed for the metagene plots.
#It uses the TAIR10_nuclear gene bed file and DM and DE gene lists generated previously
#Author: Lea Faivre
#Date: 231025

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(xlsx)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(rio)

# Data import -------------------------------------------------------------

K4Genes <- read_excel("Data/K4genes_AnyCond_CDS.xlsx")

K27Genes <- read_excel("Data/K27genes_AnyCond_CDS.xlsx")

TAIR10 <- import("Results/bedFiles/TAIR10_nuclear.bed", format = "BED")
start(TAIR10) <- start(TAIR10) - 1
#When imported into GRanges, starts positions are shifted by +1, this corrects it.

DMGenes <- import_list("Data/DMGenes_FC.xlsx", setclass = "tbl")

DEGenes <- import_list("Data/DEGenes.xlsx", setclass = "tbl")

# K4/K27 Targets ------------------------------------------------------------------------------

K4.bed <- TAIR10[TAIR10$name %in% K4Genes$Gene]
write.table(K4.bed, file = "Results/bedFiles/K4Targets.bed", 
            quote = F, sep = "\t", row.names = F, col.names = F) 

K27.bed <- TAIR10[TAIR10$name %in% K27Genes$Gene]
write.table(K27.bed, file = "Results/bedFiles/K27Targets.bed", 
            quote = F, sep = "\t", row.names = F, col.names = F) 

# DM ---------------------------------------------------------------------------------------
DM.bed <- lapply(DMGenes, function(x){
  TAIR10[TAIR10$name %in% x$Gene]
})

lapply(names(DM.bed), function(x){
  write.table(DM.bed[[x]], 
              file = paste("Results/bedFiles/DM_", x, ".bed", sep = ""), 
              quote = F, sep = "\t", row.names = F, col.names = F) 
})

K4_NonDM <- K4.bed[!(K4.bed$name %in% DMGenes$K4_3h_Gain$Gene) & 
                   !(K4.bed$name %in% DMGenes$K4_3d_Gain$Gene) &
                   !(K4.bed$name %in% DMGenes$K4_3h_Loss$Gene) &
                   !(K4.bed$name %in% DMGenes$K4_3d_Loss$Gene)]

write.table(K4_NonDM, file = "Results/bedFiles/K4_NonDM.bed", 
            quote = F, sep = "\t", row.names = F, col.names = F)

K27_NonDM <- K27.bed[!(K27.bed$name %in% DMGenes$K27_3h_Gain$Gene) & 
                     !(K27.bed$name %in% DMGenes$K27_3d_Gain$Gene) &
                     !(K27.bed$name %in% DMGenes$K27_3h_Loss$Gene) &
                     !(K27.bed$name %in% DMGenes$K27_3d_Loss$Gene)]

write.table(K27_NonDM, file = "Results/bedFiles/K27_NonDM.bed", 
            quote = F, sep = "\t", row.names = F, col.names = F)

# DE Genes (all) ------------------------------------------------------------------------------
DE <- list("h_UP" = DEGenes$Nvs3h %>%
             filter(log2FoldChange > 1),
           "h_DOWN" = DEGenes$Nvs3h %>%
             filter(log2FoldChange < -1),
           "d_UP" = DEGenes$Nvs3d %>%
             filter(log2FoldChange > 1),
           "d_DOWN" = DEGenes$Nvs3d %>%
             filter(log2FoldChange < -1))

DE.bed <- lapply(DE, function(x){
  TAIR10[TAIR10$name %in% x$Gene]
})

lapply(names(DE.bed), function(x){
  write.table(DE.bed[[x]], 
              file = paste("Results/bedFiles/DE_", x, ".bed", sep = ""), 
              quote = F, sep = "\t", row.names = F, col.names = F) 
})

NonDE <- TAIR10[!(TAIR10$name %in% DE$h_UP$Gene) &
                !(TAIR10$name %in% DE$h_DOWN$Gene) &
                !(TAIR10$name %in% DE$d_UP$Gene) &
                !(TAIR10$name %in% DE$d_DOWN$Gene)]

write.table(NonDE, file = "Results/bedFiles/NonDE.bed", 
            quote = F, sep = "\t", row.names = F, col.names = F)
