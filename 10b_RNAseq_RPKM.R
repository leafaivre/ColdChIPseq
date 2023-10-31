#This code generates the RPKM file from the featureCounts output.
#Author: Lea Faivre
#Date: 231031

# Libraries -----------------------------------------------------------------------------------
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(xlsx)
library(edgeR)

# Data import ---------------------------------------------------------------------------------
data <- read.table("InputFiles/withTE.txt", header = T, row.names = 1)

GeneLength <- data$Length

data <- data %>%
  select(-Chr, -Start, -End, -Strand, -Length)

colnames(data) <- c("N_1", "h_1", "d_1", 
                    "N_2", "h_2", "d_2", "N_3", "h_3", "d_3")

RPKM <- rpkm(data, GeneLength)

RPKM_avg <- RPKM %>%
  as.data.frame() %>%
  mutate(N = (N_1 + N_2 + N_3)/3,
         h = (h_1 + h_2 + h_3)/3,
         d = (d_1 + d_2 + d_3)/3) %>%
  select(N, h, d)

write.xlsx(RPKM_avg, "Data/RPKM_RNAseq.xlsx")
