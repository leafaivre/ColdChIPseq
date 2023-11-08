#This code uses the featureCounts file to generate fold changes and gene lists
#It uses the genes lists generated in 06a_GeneratingGenesInPeaksMasterLists.R
#Author: Lea Faivre
#Date: 231017

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(xlsx)
library(dplyr)
library(edgeR)

# Data import -------------------------------------------------------------

K4Genes <- read_excel("Data/K4genes_AnyCond_CDS.xlsx")

K27Genes <- read_excel("Data/K27genes_AnyCond_CDS.xlsx")

data <- read.table("InputFiles/FC_all_noprom.txt", header = T, row.names = 1) 

GeneLength <- data$Length

data <- data %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length)

K4 <- data %>%
  dplyr::select("N_1" = filter.A3_USR_q10.bam, "h_1" = filter.A7_USR_q10.bam,
                "d_1" = filter.A11_USR_q10.bam, "N_2" = filter.A15_USR_q10.bam,
                "h_2" = filter.A19_USR_q10.bam, "d_2" = filter.A23_USR_q10.bam)

K27 <- data %>%
  dplyr::select("N_1" = filter.A4_USR_q10.bam, "h_1" = filter.A8_USR_q10.bam,
                "d_1" = filter.A12_USR_q10.bam, "N_2" = filter.A16_USR_q10.bam,
                "h_2" = filter.A20_USR_q10.bam, "d_2" = filter.A24_USR_q10.bam)

RPKM_K4 <- rpkm(K4, GeneLength)

RPKM_avg_K4 <- RPKM_K4 %>%
  as.data.frame() %>%
  mutate(N = (N_1 + N_2)/2,
         h = (h_1 + h_2)/2,
         d = (d_1 + d_2)/2) %>%
  select(N, h, d)

#write.xlsx(RPKM_avg_K4, "Data/RPKM_K4.xlsx")

RPKM_K27 <- rpkm(K27, GeneLength)

RPKM_avg_K27 <- RPKM_K27 %>%
  as.data.frame() %>%
  mutate(N = (N_1 + N_2)/2,
         h = (h_1 + h_2)/2,
         d = (d_1 + d_2)/2) %>%
  select(N, h, d)

#write.xlsx(RPKM_avg_K27, "Data/RPKM_K27.xlsx")
