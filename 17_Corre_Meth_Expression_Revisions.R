#This code generates correlation plots between methylation levels and expression levels, as requested
#in the review process.
#It uses the RPKM files generated in 06_featureCountsRPKM.R and 10b_RNAseq_RPKM.R
#Author: Lea Faivre
#Date: 240319

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(readxl)
library(xlsx)
library(dplyr)

# Data import -------------------------------------------------------------
RPKM_Exp <- read_excel("Data/RPKM_RNAseq.xlsx")

RPKM_K27 <- read_excel("Data/RPKM_K27.xlsx")

RPKM_K4 <- read_excel("Data/RPKM_K4.xlsx")


# Correlation plots K4 ---------------------------------------------------------------------------

Cat_K4_N <- RPKM_K4 %>%
  mutate(GroupN = case_when(log10(N) <= 1 ~ "None",
                            log10(N) <= 1.5 ~ "Low",
                            log10(N) <= 2 ~ "Middle",
                            log10(N) > 2 ~ "High"),
         Grouph = case_when(log10(h) <= 1 ~ "None",
                            log10(h) <= 1.5 ~ "Low",
                            log10(h) <= 2 ~ "Middle",
                            log10(h) > 2 ~ "High"),
         Groupd = case_when(log10(d) <= 1 ~ "None",
                            log10(d) <= 1.5 ~ "Low",
                            log10(d) <= 2 ~ "Middle",
                            log10(d) > 2 ~ "High")
    
  ) %>%
  select(GroupN, Grouph, Groupd, Gene)

Exp_K4 <- inner_join(Cat_K4_N, RPKM_Exp, by = "Gene") %>%
  mutate(GroupN = fct_relevel(GroupN, "None", "Low", "Middle"),
         Grouph = fct_relevel(Grouph, "None", "Low", "Middle"),
         Groupd = fct_relevel(Groupd, "None", "Low", "Middle"))

ggplot(Exp_K4, aes(x = GroupN, y = N, color = GroupN, fill = GroupN)) +
  geom_boxplot(lwd = 0.6, alpha = 0.2) +
  scale_color_manual(values = c("black", "#8ca23d", "#717f35", "#354a12")) +
  scale_fill_manual(values = c("black", "#8ca23d", "#717f35", "#354a12")) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Gene expression (log10(RPKM))", x = "") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
ggsave("Plots/SuppFig5/K4_N.tiff", units = "in", width = 3, height = 4, dpi = 300, compression = 'lzw')

ggplot(Exp_K4, aes(x = Grouph, y = h, color = Grouph, fill = Grouph)) +
  geom_boxplot(lwd = 0.6, alpha = 0.2) +
  scale_color_manual(values = c("black", "#8ca23d", "#717f35", "#354a12")) +
  scale_fill_manual(values = c("black", "#8ca23d", "#717f35", "#354a12")) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Gene expression (log10(RPKM))", x = "") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
ggsave("Plots/SuppFig5/K4_h.tiff", units = "in", width = 3, height = 4, dpi = 300, compression = 'lzw')

ggplot(Exp_K4, aes(x = Groupd, y = d, color = Groupd, fill = Groupd)) +
  geom_boxplot(lwd = 0.6, alpha = 0.2) +
  scale_color_manual(values = c("black", "#8ca23d", "#717f35", "#354a12")) +
  scale_fill_manual(values = c("black", "#8ca23d", "#717f35", "#354a12")) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Gene expression (log10(RPKM))", x = "") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
ggsave("Plots/SuppFig5/K4_d.tiff", units = "in", width = 3, height = 4, dpi = 300, compression = 'lzw')


# Correlation plots K27 -----------------------------------------------------------------------
Cat_K27_N <- RPKM_K27 %>%
  mutate(GroupN = case_when(log10(N) <= 1 ~ "None",
                            log10(N) <= 1.5 ~ "Low",
                            log10(N) <= 2 ~ "Middle",
                            log10(N) > 2 ~ "High"),
         Grouph = case_when(log10(h) <= 1 ~ "None",
                            log10(h) <= 1.5 ~ "Low",
                            log10(h) <= 2 ~ "Middle",
                            log10(h) > 2 ~ "High"),
         Groupd = case_when(log10(d) <= 1 ~ "None",
                            log10(d) <= 1.5 ~ "Low",
                            log10(d) <= 2 ~ "Middle",
                            log10(d) > 2 ~ "High")
         
  ) %>%
  select(GroupN, Grouph, Groupd, Gene)

Exp_K27 <- inner_join(Cat_K27_N, RPKM_Exp, by = "Gene") %>%
  mutate(GroupN = fct_relevel(GroupN, "None", "Low", "Middle"),
         Grouph = fct_relevel(Grouph, "None", "Low", "Middle"),
         Groupd = fct_relevel(Groupd, "None", "Low", "Middle"))

ggplot(Exp_K27, aes(x = GroupN, y = N, color = GroupN, fill = GroupN)) +
  geom_boxplot(lwd = 0.6, alpha = 0.2) +
  scale_color_manual(values = c("black", "#bb1c02", "#8f1403", "#5b0303")) +
  scale_fill_manual(values = c("black", "#bb1c02", "#8f1403", "#5b0303")) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Gene expression (log10(RPKM))", x = "") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
ggsave("Plots/SuppFig5/K27_N.tiff", units = "in", width = 3, height = 4, dpi = 300, compression = 'lzw')

ggplot(Exp_K27, aes(x = Grouph, y = h, color = Grouph, fill = Grouph)) +
  geom_boxplot(lwd = 0.6, alpha = 0.2) +
  scale_color_manual(values = c("black", "#bb1c02", "#8f1403", "#5b0303")) +
  scale_fill_manual(values = c("black", "#bb1c02", "#8f1403", "#5b0303")) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Gene expression (log10(RPKM))", x = "") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
ggsave("Plots/SuppFig5/K27_h.tiff", units = "in", width = 3, height = 4, dpi = 300, compression = 'lzw')

ggplot(Exp_K27, aes(x = Groupd, y = d, color = Groupd, fill = Groupd)) +
  geom_boxplot(lwd = 0.6, alpha = 0.2) +
  scale_color_manual(values = c("black", "#bb1c02", "#8f1403", "#5b0303")) +
  scale_fill_manual(values = c("black", "#bb1c02", "#8f1403", "#5b0303")) +
  scale_y_log10() +
  theme_minimal() +
  labs(y = "Gene expression (log10(RPKM))", x = "") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")
ggsave("Plots/SuppFig5/K27_d.tiff", units = "in", width = 3, height = 4, dpi = 300, compression = 'lzw')

