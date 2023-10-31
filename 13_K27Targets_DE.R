#This code examines the difference between DE K27 targets which are DM or not.
#It uses the DE gene lists generated in 10_RNAseqAnalysis and 06_featureCountsGeneLists
#Author: Lea Faivre
#Date: 231031

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
library(ggpubr)

# Data import ---------------------------------------------------------------------------------
DE_FC <- list(
  h = read_excel("Data/RNAseq_Nvs3h.xlsx") %>%
    select(Gene, log2FoldChange) %>%
    dplyr::rename(DE_FC_3h = log2FoldChange),
  d = read_excel("Data/RNAseq_Nvs3d.xlsx") %>%
    select(Gene, log2FoldChange) %>%
    dplyr::rename(DE_FC_3d = log2FoldChange))

DMGenes <- import_list("Data/DMGenes_FC.xlsx", setclass = "tbl") 

DMGenes <- lapply(DMGenes, function(x){
  x %>%
    select(Gene, log2FCNvs3h, log2FCNvs3d, Name)
})

DEGenes <- import_list("Data/DEGenes.xlsx", setclass = "tbl")

DE <- list("h_UP" = DEGenes$Nvs3h %>%
             filter(log2FoldChange > 1),
           "h_DOWN" = DEGenes$Nvs3h %>%
             filter(log2FoldChange < -1),
           "d_UP" = DEGenes$Nvs3d %>%
             filter(log2FoldChange > 1),
           "d_DOWN" = DEGenes$Nvs3d %>%
             filter(log2FoldChange < -1))

K27Genes <- read_excel("Data/K27genes_AnyCond_CDS.xlsx")

RPKM <- read_excel("Data/RPKM_RNAseq.xlsx")

# Venn K27 Targets ----------------------------------------------------------------------------
##All timepoints
LosingK27 <- full_join(DMGenes$K27_3h_Loss, DMGenes$K27_3d_Loss, 
                       by = c("Gene", "log2FCNvs3h", "log2FCNvs3d", "Name" ))

GainingK27 <- full_join(DMGenes$K27_3h_Gain, DMGenes$K27_3d_Gain, 
                       by = c("Gene", "log2FCNvs3h", "log2FCNvs3d", "Name" ))

DESimp <- lapply(DE, function(x){
  x %>%
    select(Gene, Name)
})

Upregulated <- full_join(DESimp$h_UP, DESimp$d_UP, by = c("Gene", "Name"))

K27Targets <- list("H3K27me3 Targets" = K27Genes$Gene, 
             "Losing H3K27me3" = LosingK27$Gene,
             "Gaining H3K27me3" = GainingK27$Gene,
             "Upregulated" = Upregulated$Gene)

ggvenn(K27Targets, fill_color = c("#749140","#cd6155", "thistle3", "skyblue"), 
              stroke_linetype = "blank", 
              show_percentage = FALSE, fill_alpha = 0.7)
#No overrepresentation of K27 in DE genes 

##3h
K27Targets <- list("H3K27me3 Targets" = K27Genes$Gene, 
                   "Losing H3K27me3" = DMGenes$K27_3h_Loss$Gene,
                   "Upregulated" = DE$h_UP$Gene)

ggvenn(K27Targets, fill_color = c("#749140","#cd6155", "thistle3"), 
       stroke_linetype = "blank", 
       show_percentage = FALSE, fill_alpha = 0.7)
#No overrepresentation of K27 in DE genes 

##3d
K27Targets <- list("H3K27me3 Targets" = K27Genes$Gene, 
                   "Losing H3K27me3" = DMGenes$K27_3d_Loss$Gene,
                   "Upregulated" = DE$d_UP$Gene)

ggvenn(K27Targets, fill_color = c("#749140","#cd6155", "thistle3"), 
       stroke_linetype = "blank", 
       show_percentage = FALSE, fill_alpha = 0.7)
#Underrepresentation of K27 in DE genes, confirms Vyse 2020.


# Pie chart -----------------------------------------------------------------------------------
K27Targets <- list("H3K27me3 Targets" = K27Genes$Gene,
                   "Upregulated" = Upregulated$Gene)

ggvenn(K27Targets, fill_color = c("#cd6155", "thistle3"), 
       stroke_linetype = "blank", 
       show_percentage = FALSE, fill_alpha = 0.7)
#ggsave("Plots/K27Targets_UpGenes.tiff", units = "in", width = 3, height = 5, dpi = 300, compression = 'lzw')



Pie <- data.frame(Cat = c("Losing H3K27me3", "Gaining H3K27me3", "Non differentially methylated"),
                  Nb = c(20, 42, 488)) %>%
  mutate(Percentage = Nb / sum(Nb),
         Cat = as.factor(Cat),
         Label = paste(round(Percentage*100, digits = 1), "%", sep = ""))

ggplot(Pie, aes(x = "", y = Percentage, fill = Cat)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = c("darkgoldenrod","#5a8e41", "#D1AAB8")) +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "none")
#ggsave("Plots/K27_DE_Pie.tiff", units = "in", width = 3, height = 3, dpi = 300, compression = 'lzw')


# LogFC comparison ----------------------------------------------------------------------------

All <- list("Not differentially methylated" = inner_join(K27Genes, Upregulated, by = "Gene") %>%
              anti_join(GainingK27, by = "Gene") %>%
              anti_join(LosingK27, by = "Gene"), 
            "Gaining H3K27me3" = inner_join(GainingK27, Upregulated, by = c("Gene", "Name")) %>%
              select(Gene, Name), 
            "Losing H3K27me3" = inner_join(LosingK27, Upregulated, by = c("Gene", "Name")) %>%
              select(Gene, Name))

FC_Cat <- plyr::ldply(All, .id = "Category") %>%
  inner_join(DE_FC$h, by = "Gene") %>%
  inner_join(DE_FC$d, by = "Gene") 

Plot <- FC_Cat %>%
  pivot_longer(cols = c("DE_FC_3d", "DE_FC_3h"), names_to = "Timepoint", values_to = "FC") %>%
  mutate(Timepoint = recode(Timepoint, DE_FC_3h = "3h", DE_FC_3d = "3d"),
         Timepoint = as.factor(Timepoint),
         Timepoint = forcats::fct_relevel(Timepoint, "3h", "3d"))

ggplot(Plot, aes(x = Timepoint, y = FC, color = Category, fill = Category)) +
  geom_boxplot(lwd = 0.6, alpha = 0.2) +
  ylim(0, 10) +
  scale_color_manual(values = c("#945068", "goldenrod4","#365527")) +
  scale_fill_manual(values = c("#D1AAB8", "#b8860b","#5a8e41")) +
  labs(y = "Gene expression (log2 fold change)", x = "") +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
ggsave("Plots/K27_DE_LogFC.tiff", units = "in", width = 3, height = 4, dpi = 300, compression = 'lzw')

compare_means(DE_FC_3d ~ Category, data = FC_Cat)

# RPKM comparison -----------------------------------------------------------------------------
RPKM_Cat <- plyr::ldply(All, .id = "Category") %>%
  inner_join(RPKM, by = "Gene") 

Plot <- RPKM_Cat %>%
  pivot_longer(cols = c("N", "h", "d"), names_to = "Timepoint", values_to = "RPKM") %>%
  mutate(Timepoint = recode(Timepoint, h = "3h", d = "3d"),
         Timepoint = as.factor(Timepoint),
         Timepoint = forcats::fct_relevel(Timepoint, "N", "3h", "3d"))

ggplot(Plot, aes(x = Timepoint, y = log2(RPKM), color = Category, fill = Category)) +
  geom_boxplot(lwd = 0.6, alpha = 0.2) +
  scale_color_manual(values = c("#945068", "goldenrod4","#365527")) +
  scale_fill_manual(values = c("#D1AAB8", "#b8860b","#5a8e41")) +
  labs(y = "Gene expression (log2(RPKM))", x = "") +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
ggsave("Plots/K27_DE_RPKM.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

compare_means(d ~ Category, data = RPKM_Cat)
