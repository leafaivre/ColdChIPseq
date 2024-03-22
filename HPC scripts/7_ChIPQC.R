library(BiocManager)
library(ChIPQC)

setwd("/scratch/leafaivre/ColdChIP_new/results")


QC_K4 <- ChIPQC("~/ColdChIP_new/ChIPQC/K4_USR.csv")

QC_K27 <- ChIPQC("~/ColdChIP_new/ChIPQC/K27_USR.csv")

QC_H3 <- ChIPQC("~/ColdChIP_new/ChIPQC/H3_USR.csv")
save(list=ls(), file="ChIPQC_Reanalysis.Rdata")


