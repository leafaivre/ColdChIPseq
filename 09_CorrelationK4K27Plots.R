#This code plots the correlation plots of the mark levels in different conditions
#Uses the file of normalized counts generated in 06_featureCountsGeneLists
#Author: Lea Faivre
#Date: 231020


# Libraries -----------------------------------------------------------------------------------

library(readxl)
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(xlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rio)
library(plyr)


# Data import ---------------------------------------------------------------------------------

Data <- import_list("Data/FoldChanges_allGenes.xlsx", setclass = "tbl")


# Graphs --------------------------------------------------------------------------------------

ggplot(Data$K4, aes(x = N, y = h)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  xlim(0, 8000) +
  ylim(0, 8000) +
  theme_minimal() +
  geom_segment(data = Data$K4, aes(x = 0, xend = 8000, y = 0, yend = 8000),
               color = "darkgreen", linewidth = 0.8) +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "Norm. counts 3h", x = "Norm. counts N")
ggsave("Plots/K4_Nvs3h.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

ggplot(Data$K4, aes(x = N, y = d)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  xlim(0, 8000) +
  ylim(0, 8000) +
  theme_minimal() +
  geom_segment(data = Data$K4, aes(x = 0, xend = 8000, y = 0, yend = 8000), 
               color = "darkgreen", linewidth = 0.8) +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "Norm. counts 3d", x = "Norm. counts N")
ggsave("Plots/K4_Nvs3d.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

ggplot(Data$K27, aes(x = N, y = h)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  xlim(0, 10000) +
  ylim(0, 10000) +
  theme_minimal() +
  geom_segment(data = Data$K4, aes(x = 0, xend = 10000, y = 0, yend = 10000),
               color = "darkred", linewidth = 0.8) +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "Norm. counts 3h", x = "Norm. counts N")
ggsave("Plots/K27_Nvs3h.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

ggplot(Data$K27, aes(x = N, y = d)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  xlim(0, 10000) +
  ylim(0, 10000) +
  theme_minimal() +
  geom_segment(data = Data$K4, aes(x = 0, xend = 10000, y = 0, yend = 10000), 
               color = "darkred", linewidth = 0.8) +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "Norm. counts 3d", x = "Norm. counts N")
ggsave("Plots/K27_Nvs3d.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')


# Cross correlation ---------------------------------------------------------------------------
OnlyRelevant <- lapply(Data, function(x){
  x %>%
    select(Gene, log2FCNvs3h, log2FCNvs3d)
})

K4K27 <- inner_join(OnlyRelevant$K4, OnlyRelevant$K27, by = "Gene", suffix = c(".K4", ".K27"))

ggplot(K4K27, aes(x = log2FCNvs3d.K4, y = log2FCNvs3d.K27)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  theme_minimal() +
  geom_smooth(method = lm) +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "K27", x = "K4")

##No correlation between H3K4me3 and H3K27me3 changes