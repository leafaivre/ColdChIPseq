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
library(plyranges)

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
    filter(abs(log2FC) > 0.5) %>%
<<<<<<< HEAD
    select(Chrom, Start, End, Event, log2FC, padj)
=======
    select(Chrom, Start, End, Length, Event, log2FC, padj)
>>>>>>> 9844bf5ad70849fa555bd91841a7a768896f820e
})

## Converting a GRange

GR <- lapply(diffRep, function(x){
  makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
})

## Filtering on the background: the DM should belong to a peak 
InPeak <- list(
  K4_3h = GR$K4_3h[overlapsAny(GR$K4_3h, Peaks_K4_3h, minoverlap = 50)],
  K4_3d = GR$K4_3d[overlapsAny(GR$K4_3d, Peaks_K4_3d, minoverlap = 50)],
  K27_3h = GR$K27_3h[overlapsAny(GR$K27_3h, Peaks_K27_3h, minoverlap = 200)],
  K27_3d = GR$K27_3d[overlapsAny(GR$K27_3d, Peaks_K27_3d, minoverlap = 200)])

## Dividing on UP/DOWN

UDRegions <- list(
  K4hU = InPeak$K4_3h %>%
    plyranges::filter(Event == "Up"),
  K4hD = InPeak$K4_3h %>%
    plyranges::filter(Event == "Down"),
  K4dU = InPeak$K4_3d %>%
    plyranges::filter(Event == "Up"),
  K4dD = InPeak$K4_3d %>%
    plyranges::filter(Event == "Down"),
  K27hU = InPeak$K27_3h %>%
    plyranges::filter(Event == "Up"),
  K27hD = InPeak$K27_3h %>%
    plyranges::filter(Event == "Down"),
  K27dU = InPeak$K27_3d %>%
    plyranges::filter(Event == "Up"),
  K27dD = InPeak$K27_3d %>%
    plyranges::filter(Event == "Down"))



# Creating bed tracks -----------------------------------------------------

lapply(names(UDRegions), function(x){
  write.table(UDRegions[[x]], 
              file = paste("Data/", x, ".bed", sep = ""), 
              quote = F, sep = "\t", row.names = F, col.names = F) 
})
