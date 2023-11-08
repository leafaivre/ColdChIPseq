#This code examines the relationships between DM and DE genes
#It uses the DE gene lists generated in 10_RNAseqAnalysis and 06_featureCountsGeneLists
#Author: Lea Faivre
#Date: 231024


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

# Venn diagrams -------------------------------------------------------------------------------
h_UP <- list("Gaining H3K4me3" = DMGenes$K4_3h_Gain$Gene, 
             "Losing H3K27me3" = DMGenes$K27_3h_Loss$Gene,
             "Upregulated" = DE$h_UP$Gene)

h_U <- ggvenn(h_UP, fill_color = c("#749140","#cd6155", "thistle3"), 
       stroke_linetype = "blank", 
       show_percentage = FALSE, fill_alpha = 0.7)

h_DOWN <- list("Losing H3K4me3" = DMGenes$K4_3h_Loss$Gene, 
               "Gaining H3K27me3" = DMGenes$K27_3h_Gain$Gene,
               "Downregulated" = DE$h_DOWN$Gene)

h_D <- ggvenn(h_DOWN, fill_color = c("#749140","#cd6155", "powderblue"), 
              stroke_linetype = "blank", 
              show_percentage = FALSE, fill_alpha = 0.7)

d_UP <- list("Gaining H3K4me3" = DMGenes$K4_3d_Gain$Gene, 
             "Losing H3K27me3" = DMGenes$K27_3d_Loss$Gene,
             "Upregulated" = DE$d_UP$Gene)

d_U <- ggvenn(d_UP, fill_color = c("#749140","#cd6155", "thistle3"), 
              stroke_linetype = "blank", 
              show_percentage = FALSE, fill_alpha = 0.7)

d_DOWN <- list("Losing H3K4me3" = DMGenes$K4_3d_Loss$Gene, 
               "Gaining H3K27me3" = DMGenes$K27_3d_Gain$Gene,
               "Downregulated" = DE$d_DOWN$Gene)

d_D <- ggvenn(d_DOWN, fill_color = c("#749140","#cd6155", "powderblue"), 
       stroke_linetype = "blank", 
       show_percentage = FALSE, fill_alpha = 0.7)


grid.newpage()
cowplot::plot_grid(h_U, h_D, d_U, d_D, 
          labels = c("Up 3h", "Down 3h", "Up 3d", "Down 3d"), 
          label_x = 0.5, hjust = 0.5)
ggsave("Plots/DEDM_Overlaps.tiff", units = "in", width = 6, height = 6, dpi = 300, compression = 'lzw')


# Correlations all genes --------------------------------------------------------------------------------

FullCorr <- full_join(DE_FC$h, DE_FC$d, by = "Gene") %>%
  full_join(DM_FC_form$K4, by = "Gene") %>%
  full_join(DM_FC_form$K27, by = c("Gene"))


ggplot(FullCorr, aes(x = K4_FC_3h, y =  DE_FC_3h)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  geom_smooth(method = lm)  +
  ylim(-10, 10) +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "log2FC(3h - N) in RNAseq", x = "log2FC(3h - N) in H3K4me3 ChIPseq")
#ggsave("Plots/Correlations/All_K4_RNA_3h.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

ggplot(FullCorr, aes(x = K4_FC_3d, y =  DE_FC_3d)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  geom_smooth(method = lm)  +
  ylim(-10, 10) +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "log2FC(3d - N) in RNAseq", x = "log2FC(3d - N) in H3K4me3 ChIPseq")
#ggsave("Plots/Correlations/All_K4_RNA_3d.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

ggplot(FullCorr, aes(x = K27_FC_3h, y =  DE_FC_3h)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  geom_smooth(method = lm)  +
  #ylim(-10, 10) +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "log2FC(3h - N) in RNAseq", x = "log2FC(3h - N) in H3K27me3 ChIPseq")
#ggsave("Plots/Correlations/All_K27_RNA_3h.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

ggplot(FullCorr, aes(x = K27_FC_3d, y =  DE_FC_3d)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  geom_smooth(method = lm)  +
  #ylim(-10, 10) +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "log2FC(3d - N) in RNAseq", x = "log2FC(3d - N) in H3K27me3 ChIPseq")
#ggsave("Plots/Correlations/All_K27_RNA_3d.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

cor.test(FullCorr$DE_FC_3h, FullCorr$K4_FC_3h, method = "spearman")
cor.test(FullCorr$DE_FC_3h, FullCorr$K27_FC_3h, method = "spearman")
cor.test(FullCorr$DE_FC_3d, FullCorr$K4_FC_3d, method = "spearman")
cor.test(FullCorr$DE_FC_3d, FullCorr$K27_FC_3d, method = "spearman")

# Correlation DE genes ------------------------------------------------------------------------

hCorr <- FullCorr %>%
  filter(padj_3h < 0.05)

dCorr <- FullCorr %>%
  filter(padj_3d < 0.05)

ggplot(hCorr, aes(x = K4_FC_3h, y =  DE_FC_3h)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  geom_smooth(method = lm)  +
  ylim(-10, 10) +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "log2FC(3h - N) in RNAseq", x = "log2FC(3h - N) in H3K4me3 ChIPseq")
#ggsave("Plots/Correlations/DE_K4_RNA_3h.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

ggplot(dCorr, aes(x = K4_FC_3d, y =  DE_FC_3d)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  geom_smooth(method = lm)  +
  ylim(-10, 10) +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "log2FC(3d - N) in RNAseq", x = "log2FC(3d - N) in H3K4me3 ChIPseq")
#ggsave("Plots/Correlations/DE_K4_RNA_3d.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

ggplot(hCorr, aes(x = K27_FC_3h, y =  DE_FC_3h)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  geom_smooth(method = lm)  +
  ylim(-10, 10) +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "log2FC(3h - N) in RNAseq", x = "log2FC(3h - N) in H3K27me3 ChIPseq")
#ggsave("Plots/Correlations/DE_K27_RNA_3h.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

ggplot(dCorr, aes(x = K27_FC_3d, y =  DE_FC_3d)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_vline(xintercept = 0, color = "darkgray") +
  geom_point(alpha = 0.2) +
  geom_smooth(method = lm)  +
  ylim(-10, 10) +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "log2FC(3d - N) in RNAseq", x = "log2FC(3d - N) in H3K27me3 ChIPseq")
#ggsave("Plots/Correlations/DE_K27_RNA_3d.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')

cor.test(hCorr$DE_FC_3h, hCorr$K4_FC_3h, method = "spearman")
cor.test(hCorr$DE_FC_3h, hCorr$K27_FC_3h, method = "spearman")
cor.test(dCorr$DE_FC_3d, dCorr$K4_FC_3d, method = "spearman")
cor.test(dCorr$DE_FC_3d, dCorr$K27_FC_3d, method = "spearman")


# Cross-timepoints ----------------------------------------------------------------------------

cor.test(hCorr$DE_FC_3h, hCorr$K4_FC_3d, method = "spearman")
cor.test(hCorr$DE_FC_3h, hCorr$K27_FC_3d, method = "spearman")
cor.test(dCorr$DE_FC_3d, dCorr$K4_FC_3h, method = "spearman")
cor.test(dCorr$DE_FC_3d, dCorr$K27_FC_3h, method = "spearman")
