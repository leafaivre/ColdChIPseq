#This code uses the featureCounts file to generate fold changes and gene lists
#It uses the genes lists generated in 06a_GeneratingGenesInPeaksMasterLists.R
#Author: Lea Faivre
#Date: 231017

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

TAIR10genes <- read_excel("InputFiles/TAIR10genes.xlsx", 
                          sheet = "Name", col_types = c("text", "text"))

K4Genes <- read_excel("Data/K4genes_AnyCond_CDS.xlsx")

K27Genes <- read_excel("Data/K27genes_AnyCond_CDS.xlsx")

data <- read.table("InputFiles/FC_all_noprom.txt", header = T, row.names = 1) %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length)

K4 <- data %>%
  dplyr::select("N_1" = filter.A3_USR_q10.bam, "h_1" = filter.A7_USR_q10.bam,
                "d_1" = filter.A11_USR_q10.bam, "N_2" = filter.A15_USR_q10.bam,
                "h_2" = filter.A19_USR_q10.bam, "d_2" = filter.A23_USR_q10.bam)

K27 <- data %>%
  dplyr::select("N_1" = filter.A4_USR_q10.bam, "h_1" = filter.A8_USR_q10.bam,
                "d_1" = filter.A12_USR_q10.bam, "N_2" = filter.A16_USR_q10.bam,
                "h_2" = filter.A20_USR_q10.bam, "d_2" = filter.A24_USR_q10.bam)


meta <- read.table("InputFiles/metadata.txt", header = T, row.names = 1)
meta <- meta %>% 
  mutate(Rep = as.factor(Rep))



data <- list(K4 = K4, K27 = K27)

dds <- lapply(data, function(x){
  DESeqDataSetFromMatrix(countData = x, 
                         colData = meta, 
                         design = ~ Rep + Condition)
})


# Normalization -----------------------------------------------------------

dds <- lapply(dds, estimateSizeFactors)
normalized_counts <- lapply(dds, function(x){
  counts(x, normalized = TRUE)
})

NormCounts <- lapply(normalized_counts, function(x){
  x %>%
    data.frame() %>%
    rownames_to_column(var = "Gene") %>%
    as_tibble()
})


# Getting the log2FC ------------------------------------------------------

NormCountsAvg <- lapply(NormCounts, function(x){
  x %>%
    mutate(N = (N_1 + N_2)/2,
           h = (h_1 + h_2)/2,
           d = (d_1 + d_2)/2,
           log2FCNvs3h = log2(h/N),
           log2FCNvs3d = log2(d/N)) 
})

NormCountsAvg <- lapply(NormCountsAvg, function(x){
  x %>%
    left_join(TAIR10genes, by  = "Gene") 
})

# Generating the files with all fold changes ------------------------------

#write.xlsx(as.data.frame(NormCountsAvg$K4),"Data/FC_featureCount_K4.xlsx", 
           #row.names = FALSE)

#write.xlsx(as.data.frame(NormCountsAvg$K27),"Data/FC_featureCount_K27.xlsx", 
           #row.names = FALSE)


# Generating the gene lists -----------------------------------------------

## Keeping only the genes in peaks
InPeaks <- list(
  K4 = inner_join(NormCountsAvg$K4, K4Genes, by = "Gene"),
  K27 = inner_join(NormCountsAvg$K27, K27Genes, by = "Gene"))

## Keeping only the significant
InPeaksFC <- list(
  K4_3h_Gain = InPeaks$K4 %>%
    filter(log2FCNvs3h > 0.5) %>%
    filter(N_1 < h_1 & N_2 < h_2),
  
  K4_3h_Loss = InPeaks$K4 %>%
    filter(log2FCNvs3h < -0.5) %>%
    filter(N_1 > h_1 & N_2 > h_2),
  
  K4_3d_Gain = InPeaks$K4 %>%
    filter(log2FCNvs3d > 0.5) %>%
    filter(N_1 < d_1 & N_2 < d_2),
  
  K4_3d_Loss = InPeaks$K4 %>%
    filter(log2FCNvs3d < -0.5) %>%
    filter(N_1 > d_1 & N_2 > d_2),
  
  K27_3h_Gain = InPeaks$K27 %>%
    filter(log2FCNvs3h > 0.5) %>%
    filter(N_1 < h_1 & N_2 < h_2),
  
  K27_3h_Loss = InPeaks$K27 %>%
    filter(log2FCNvs3h < -0.5) %>%
    filter(N_1 > h_1 & N_2 > h_2),
  
  K27_3d_Gain = InPeaks$K27 %>%
    filter(log2FCNvs3d > 0.5) %>%
    filter(N_1 < d_1 & N_2 < d_2),
  
  K27_3d_Loss = InPeaks$K27 %>%
    filter(log2FCNvs3d < -0.5) %>%
    filter(N_1 > d_1 & N_2 > d_2))


# Writing the files -------------------------------------------------------

for (i in names(InPeaksFC)) {
    write.xlsx(x = as.data.frame(InPeaksFC[[i]]),
               file = "Data/DMGenes_FC.xlsx",
               sheetName = paste(i, sep = ""),
               row.names = FALSE,
               append = T)
  }

