#This code characterizes the significant diffReps regions
#Uses the tracks generated in 01_GeneratingFilteredLists
#Author: Lea Faivre
#Date: 231011


# Libraries ---------------------------------------------------------------

library(plyr)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(readxl)
library(plyranges)
library(ggpubr)
library(rcompanion)
library(FSA)
library(ChIPseeker)

Araport <- TxDb.Athaliana.BioMart.plantsmart28 

# Data import -------------------------------------------------------------
extraCols_FC <- c(Length = "numeric", Strand = "character", Event = "character", log2FC = "numeric", padj = "numeric")

DMR_K4 <- list(
  "3h Gain" = rtracklayer::import("Data/K4hU.bed", format = "BED", 
                                  extraCols = extraCols_FC),
  "3h Loss" = rtracklayer::import("Data/K4hD.bed", format = "BED", 
                                  extraCols = extraCols_FC),
  "3d Gain" = rtracklayer::import("Data/K4dU.bed", format = "BED",
                                  extraCols = extraCols_FC),
  "3d Loss" = rtracklayer::import("Data/K4dD.bed", format = "BED", 
                                  extraCols = extraCols_FC))
DMR_K27 <- list(
  "3h Gain" = rtracklayer::import("Data/K27hU.bed", format = "BED", 
                                  extraCols = extraCols_FC),
  "3h Loss" = rtracklayer::import("Data/K27hD.bed", format = "BED", 
                                  extraCols = extraCols_FC),
  "3d Gain" = rtracklayer::import("Data/K27dU.bed", format = "BED", 
                                  extraCols = extraCols_FC),
  "3d Loss" = rtracklayer::import("Data/K27dD.bed", format = "BED", 
                                  extraCols = extraCols_FC))

# Plotting the length -----------------------------------------------------
DMR_K4_length <- lapply(DMR_K4, width) %>%
  ldply(data.frame, .id = "Sample") %>%
  dplyr::rename(Length = 'X..i..') %>%
  dplyr::mutate(Sample = as.factor(Sample),
                Sample = forcats::fct_relevel(Sample, "3h Gain", "3h Loss", 
                                              "3d Gain", "3d Loss"))

ggplot(DMR_K4_length, aes(x = Sample, y = Length)) +
  geom_violin(aes(x = Sample, y = Length, fill = Sample), scale = "area", 
              alpha = 0.7, linetype = "blank") +
  geom_boxplot(aes(x = Sample, y = Length), width = 0.1, fill = "white", 
               outlier.shape = NA) +
  ylim(50, 300) +
  scale_fill_manual(values = c("#749140","#cd6155","#196F3D","#5f0000")) +
  labs(y = "Length (bp)", x = "") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/K4DMRLength.tiff", units = "in", width = 3, height = 4, dpi = 300,
       #compression = 'lzw')

DMR_K27_length <- lapply(DMR_K27, width) %>%
  ldply(data.frame, .id = "Sample") %>%
  dplyr::rename(Length = 'X..i..') %>%
  dplyr::mutate(Sample = as.factor(Sample),
                Sample = forcats::fct_relevel(Sample, "3h Gain", "3h Loss", 
                                              "3d Gain", "3d Loss"))

ggplot(DMR_K27_length, aes(x = Sample, y = Length)) +
  geom_violin(aes(x = Sample, y = Length, fill = Sample), scale = "area", 
              alpha = 0.7, linetype = "blank") +
  geom_boxplot(aes(x = Sample, y = Length), width = 0.1, fill = "white", 
               outlier.shape = NA) +
  ylim(200, 1000) +
  scale_fill_manual(values = c("#749140","#cd6155","#196F3D","#5f0000")) +
  labs(y = "Length (bp)", x = "") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/K27DMRLength.tiff", units = "in", width = 3, height = 4, dpi = 300, 
       #compression = 'lzw')

# Plotting the fold change ------------------------------------------------
DMR_K4_FC <- lapply(DMR_K4, function(x){
  x$log2FC }) %>%
  ldply(data.frame, .id = "Sample") %>%
  dplyr::rename(log2FC = 'X..i..') %>%
  dplyr::mutate(Sample = as.factor(Sample),
                Sample = forcats::fct_relevel(Sample, "3h Gain", "3h Loss", 
                                              "3d Gain", "3d Loss"),
                log2FC = abs(log2FC))

ggplot(DMR_K4_FC, aes(x = Sample, y = log2FC)) +
  geom_violin(aes(x = Sample, y = log2FC, fill = Sample), scale = "area", 
              alpha = 0.7, linetype = "blank") +
  geom_boxplot(aes(x = Sample, y = log2FC), width = 0.1, fill = "white", 
               outlier.shape = NA) +
  scale_fill_manual(values = c("#749140","#cd6155","#196F3D","#5f0000") ) +
  labs(y = "|log2FC|", x = "") +
  ylim(0.3, 2) +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/K4DMRFC.tiff", units = "in", width = 6, height = 4, dpi = 300, 
       #compression = 'lzw')

DMR_K27_FC <- lapply(DMR_K27, function(x){
  x$log2FC }) %>%
  ldply(data.frame, .id = "Sample") %>%
  dplyr::rename(log2FC = 'X..i..') %>%
  dplyr::mutate(Sample = as.factor(Sample),
                Sample = forcats::fct_relevel(Sample, "3h Gain", "3h Loss", 
                                              "3d Gain", "3d Loss"),
                log2FC = abs(log2FC))

ggplot(DMR_K27_FC, aes(x = Sample, y = log2FC)) +
  geom_violin(aes(x = Sample, y = log2FC, fill = Sample), scale = "area", 
              alpha = 0.7, linetype = "blank") +
  geom_boxplot(aes(x = Sample, y = log2FC), width = 0.1, fill = "white", 
               outlier.shape = NA) +
  scale_fill_manual(values = c("#749140","#cd6155","#196F3D","#5f0000") ) +
  labs(y = "|log2FC|", x = "") +
  ylim(0.3, 2) +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
#ggsave("Plots/K27DMRFC.tiff", units = "in", width = 6, height = 4, dpi = 300, 
       #compression = 'lzw')

# Annotation --------------------------------------------------------------
AnnoK4 <- lapply(DMR_K4, annotatePeak, TxDb = Araport, level = "gene",
                 genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron", 
                                               "Promoter", "Downstream", "Intergenic"),
               tssRegion = c(-500, 500), verbose = FALSE)

plotAnnoBar(AnnoK4)
#ggsave("Plots/K4DMRAnno.tiff", units = "in", width = 8, height = 4, dpi = 300, 
       #compression = 'lzw')
plotDistToTSS(AnnoK4)

AnnoK27 <- lapply(DMR_K27, annotatePeak, TxDb = Araport, level = "gene",
                 genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron", 
                                               "Promoter", "Downstream", "Intergenic"),
                 tssRegion = c(-500, 500), verbose = FALSE)

plotAnnoBar(AnnoK27)
#ggsave("Plots/K27DMRAnno.tiff", units = "in", width = 8, height = 4, dpi = 300, 
       #compression = 'lzw')
plotDistToTSS(AnnoK27)