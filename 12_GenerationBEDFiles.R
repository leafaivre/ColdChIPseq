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

TAIR10 <- rtracklayer::import("Results/bedFiles/TAIR10_nuclear.bed", format = "BED")

DMGenes <- import_list("Data/DMGenes_FC.xlsx", setclass = "tbl")

DEGenes <- import_list("Data/DEGenes.xlsx", setclass = "tbl")

# K4/K27 Targets ------------------------------------------------------------------------------

K4.bed <- TAIR10[TAIR10$name %in% K4Genes$Gene]
#export.bed(K4.bed, con = "Results/bedFiles/K4Targets.bed")

K27.bed <- TAIR10[TAIR10$name %in% K27Genes$Gene]
#export.bed(K27.bed, con = "Results/bedFiles/K27Targets.bed")

# DM ---------------------------------------------------------------------------------------
DM.bed <- lapply(DMGenes, function(x){
  TAIR10[TAIR10$name %in% x$Gene]
})

lapply(names(DM.bed), function(x){
  #export.bed(DM.bed[[x]], con = paste("Results/bedFiles/DM_", x, ".bed", sep = ""))
})

K4_NonDM <- K4.bed[!(K4.bed$name %in% DMGenes$K4_3h_Gain$Gene) & 
                   !(K4.bed$name %in% DMGenes$K4_3d_Gain$Gene) &
                   !(K4.bed$name %in% DMGenes$K4_3h_Loss$Gene) &
                   !(K4.bed$name %in% DMGenes$K4_3d_Loss$Gene)]

#export.bed(K4_NonDM, con = "Results/bedFiles/K4_NonDM.bed",)

K27_NonDM <- K27.bed[!(K27.bed$name %in% DMGenes$K27_3h_Gain$Gene) & 
                     !(K27.bed$name %in% DMGenes$K27_3d_Gain$Gene) &
                     !(K27.bed$name %in% DMGenes$K27_3h_Loss$Gene) &
                     !(K27.bed$name %in% DMGenes$K27_3d_Loss$Gene)]

#export.bed(K27_NonDM, con = "Results/bedFiles/K27_NonDM.bed")

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
  #export.bed(DE.bed[[x]], con = paste("Results/bedFiles/DE_", x, ".bed", sep = "")) 
})

NonDE <- TAIR10[!(TAIR10$name %in% DE$h_UP$Gene) &
                !(TAIR10$name %in% DE$h_DOWN$Gene) &
                !(TAIR10$name %in% DE$d_UP$Gene) &
                !(TAIR10$name %in% DE$d_DOWN$Gene)]

#export.bed(NonDE, con = "Results/bedFiles/NonDE.bed")

# DE Genes (Targets only) ---------------------------------------------------------------------

DE.K4.bed <- lapply(DE, function(x){
  K4.bed[K4.bed$name %in% x$Gene]
})

lapply(names(DE.K4.bed), function(x){
  #export.bed(DE.K4.bed[[x]], con = paste("Results/bedFiles/DE_K4", x, ".bed", sep = "")) 
})

NonDEK4 <- K4.bed[!(K4.bed$name %in% DE$h_UP$Gene) &
                  !(K4.bed$name %in% DE$h_DOWN$Gene) &
                  !(K4.bed$name %in% DE$d_UP$Gene) &
                  !(K4.bed$name %in% DE$d_DOWN$Gene)]

#export.bed(NonDEK4, con = "Results/bedFiles/NonDE_K4.bed")

DE.K27.bed <- lapply(DE, function(x){
  K27.bed[K27.bed$name %in% x$Gene]
})

lapply(names(DE.K27.bed), function(x){
  #export.bed(DE.K27.bed[[x]], con = paste("Results/bedFiles/DE_K27", x, ".bed", sep = "")) 
})

NonDEK27 <- K27.bed[!(K27.bed$name %in% DE$h_UP$Gene) &
                    !(K27.bed$name %in% DE$h_DOWN$Gene) &
                    !(K27.bed$name %in% DE$d_UP$Gene) &
                    !(K27.bed$name %in% DE$d_DOWN$Gene)]

#export.bed(NonDEK4, con = "Results/bedFiles/NonDE_K4.bed")
