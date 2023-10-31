#This code creates lists of genes in the H3K4me3 (resp. H3K27me3) peaks
# Uses the master peaks lists generated in 01a
#Author: Lea Faivre
#Date: 231017

# Libraries ---------------------------------------------------------------
library(GenomicRanges)
library(dplyr)
library(plyr)
library(ggplot2)
library(tidyr)
library(rtracklayer)
library(HelloRanges)
library(TxDb.Athaliana.BioMart.plantsmart28)

# Import and format data -------------------------------------------------------
Peaks_K4 <- rtracklayer::import("Data/K4.N.3h.3d.peaks.bed", format = "BED")

Peaks_K27 <- rtracklayer::import("Data/K27.N.3h.3d.peaks.bed", format = "BED")

NPeaks_K27 <- read.delim("InputFiles/A4A16-idr", header = FALSE)

Araport <- TxDb.Athaliana.BioMart.plantsmart28

Genes <- genes(Araport, columns = c("GENEID", "TXTYPE"))


# Getting the K27 genes -------------------------------------------------------

K27genes <- subsetByOverlaps(Genes, Peaks_K27, minoverlap = 150, ignore.strand = TRUE)

K27.genes <- as.data.frame(unlist(K27genes$GENEID))

colnames(K27.genes) <- ("Gene")

xlsx::write.xlsx(K27.genes, "Data/K27genes_AnyCond_CDS.xlsx", row.names = FALSE)


# Getting the K4 genes ----------------------------------------------------

K4genes <- subsetByOverlaps(Genes, Peaks_K4, minoverlap = 150, ignore.strand = TRUE)

K4.genes <- as.data.frame(unlist(K4genes$GENEID))

colnames(K4.genes) <- ("Gene")

xlsx::write.xlsx(K4.genes, "Data/K4genes_AnyCond_CDS.xlsx", row.names = FALSE)


# K27 Genes in N ------------------------------------------------------------------------------
colnames(NPeaks_K27) <- c("Chr", "Start", "End", "Name", "IDR", "Strand", "Signal", 
                 "pvalue", "qvalue", "localIDR", "globalIDR", 
                 "R1Start", "R1End", "R1pvalue",
                 "R2Start", "R2End", "R2pvalue")

NPeaks_K27GR <- makeGRangesFromDataFrame(NPeaks_K27, keep.extra.columns = TRUE, 
                         start.field = "Start", end.field = "End")



NK27genes <- subsetByOverlaps(Genes, NPeaks_K27GR, minoverlap = 150, ignore.strand = TRUE)

NK27.genes <- as.data.frame(unlist(NK27genes$GENEID))

colnames(NK27.genes) <- ("Gene")

xlsx::write.xlsx(NK27.genes, "Data/K27genes_N_CDS.xlsx", row.names = FALSE)