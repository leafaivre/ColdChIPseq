#This code creates master lists for peaks in different conditions based on the 
#IDR results
#Author: Lea Faivre
#Date: 231009


# Libraries ---------------------------------------------------------------
library(GenomicRanges)
library(dplyr)
library(plyr)
library(ggplot2)
library(tidyr)
library(rtracklayer)
library(HelloRanges)


# Import and format data -------------------------------------------------------
AllIDR <- list(
  N.K27 = read.delim("InputFiles/A4A16-idr", header = FALSE),
  h.K27 = read.delim("InputFiles/A8A20-idr", header = FALSE),
  d.K27 = read.delim("InputFiles/A12A24-idr", header = FALSE),
  N.K4 = read.delim("InputFiles/A3A15-idr", header = FALSE),
  h.K4 = read.delim("InputFiles/A7A19-idr", header = FALSE),
  d.K4 = read.delim("InputFiles/A11A23-idr", header = FALSE))

ChangeNames <- function(x){
  colnames(x) <- c("Chr", "Start", "End", "Name", "IDR", "Strand", "Signal", 
                       "pvalue", "qvalue", "localIDR", "globalIDR", 
                       "R1Start", "R1End", "R1pvalue",
                       "R2Start", "R2End", "R2pvalue")
  
  makeGRangesFromDataFrame(x, keep.extra.columns = TRUE, 
                           start.field = "Start",
                           end.field = "End")
}

GR_IDR <- lapply(AllIDR, ChangeNames)


# Creating the master lists -----------------------------------------------

K4_3hmaster <- GenomicRanges::union(GR_IDR$N.K4, GR_IDR$h.K4) 
export.bed(K4_3hmaster, "Data/K4.N.3h.peaks.bed")

K4_3dmaster <- GenomicRanges::union(GR_IDR$N.K4, GR_IDR$d.K4) 
export.bed(K4_3dmaster, "Data/K4.N.3d.peaks.bed")

AllK4 <- GenomicRanges::union(K4_3hmaster, K4_3dmaster)
export.bed(AllK4, "Data/K4.N.3h.3d.peaks.bed")

K27_3hmaster <- GenomicRanges::union(GR_IDR$N.K27, GR_IDR$h.K27) 
export.bed(K27_3hmaster, "Data/K27.N.3h.peaks.bed")

K27_3dmaster <- GenomicRanges::union(GR_IDR$N.K27, GR_IDR$d.K27) 
export.bed(K27_3dmaster, "Data/K27.N.3d.peaks.bed")

AllK27 <- GenomicRanges::union(K27_3hmaster, K27_3dmaster)
export.bed(AllK27, "Data/K27.N.3h.3d.peaks.bed")