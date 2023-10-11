#This code filters the diffReps output files to keep only the significantly differentially methylated regions
#Author: Lea Faivre
#Date: 231009


# Libraries ---------------------------------------------------------------

library(plyr)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(readxl)


# Data import -------------------------------------------------------------
diffRep <- list(K4_3h = read.delim("InputFiles/diff.K24.3h.50.midfrag.bed"),
                K4_3d = read.delim("InputFiles/diff.K24.3d.50.midfrag.bed"),
                K27_3h = read.delim("InputFiles/diff.K27.3h.200.midfrag.bed"),
                K27_3d = read.delim("InputFiles/diff.K27.3d.200.midfrag.bed"))

Peaks_K4_3h <- rtracklayer::import("Data/K4.N.3h.peaks.bed", format = "BED")

Peaks_K4_3d <- rtracklayer::import("Data/K4.N.3d.peaks.bed", format = "BED")

Peaks_K27_3h <- rtracklayer::import("Data/K27.N.3h.peaks.bed", format = "BED")

Peaks_K27_3d <- rtracklayer::import("Data/K27.N.3d.peaks.bed", format = "BED")

# Filtering ---------------------------------------------------------------

## Filtering is done on the p-value and on the fold change (log2FC > 2), which 
## represents a gain of at least ~50% or a loss of at least ~30%. 

diffRep <- lapply(diffRep, function(x){
  x %>%
    filter(!(Chrom %in% c("Pt", "Mt"))) %>%
    filter(padj < 0.001) %>%
    filter(abs(log2FC) > 0.5)
})

## Converting a GRange

GR <- lapply(diffRep, function(x){
  makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
})

## Filtering on the background: the DM should overlap a peak 
InPeak <- list(
  K4_3h = GR$K4_3h[overlapsAny(GR$K4_3h, Peaks_K4_3h, minoverlap = 1)],
  K4_3d = GR$K4_3d[overlapsAny(GR$K4_3d, Peaks_K4_3d, minoverlap = 1)],
  K27_3h = GR$K27_3h[overlapsAny(GR$K27_3h, Peaks_K27_3h, minoverlap = 1)],
  K27_3d = GR$K27_3d[overlapsAny(GR$K27_3d, Peaks_K27_3d, minoverlap = 1)])

## Dividing on UP/DOWN

UDRegions <- list(
  K4hU = InPeak$K4_3h %>%
)
plyranges::filter


UDRegions <- list(
  K43hU = diffRep$K4_3h %>%
    filter(Event == "Up"),
  K43hD = diffRep$K4_3h %>%
    filter(Event == "Down"),
  K43dU = diffRep$K4_3d %>%
    filter(Event == "Up"),
  K43dD = diffRep$K4_3d %>%
    filter(Event == "Down"),
  K273hU = diffRep$K27_3h %>%
    filter(Event == "Up"),
  K273hD = diffRep$K27_3h %>%
    filter(Event == "Down"),
  K273dU = diffRep$K27_3d %>%
    filter(Event == "Up"),
  K273dD = diffRep$K27_3d %>%
    filter(Event == "Down"))




