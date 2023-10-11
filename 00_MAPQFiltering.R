#This code examines the impact on the MAPQ filtering on the detection of differentially methylated genes
#Author: Lea Faivre
#Date: 231005


# Libraries and necessary files ------------------------------------------------

library(plyr)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(readxl)
library(ggvenn)

Araport <- TxDb.Athaliana.BioMart.plantsmart28 

Genes <- genes(Araport, columns = c("GENEID","TXTYPE"))

TAIR10genes <- read_excel("InputFiles/TAIR10genes.xlsx", 
                          sheet = "Name", col_types = c("text", "text"))


# Data import and formatting ----------------------------------------------

Filt <- read.delim("InputFiles/diff.K27.3d.200.midfrag.bed") %>%
  filter(!(Chrom %in% c("Pt", "Mt"))) %>%
  filter(padj < 0.001) %>%
  filter(log2FC < -0.322)  %>%
  dplyr::select(Chrom, Start, End, Length, Event, log2FC, pval, padj)

GRFilt <- makeGRangesFromDataFrame(Filt, keep.extra.columns = TRUE)

UnFilt <- read.delim("InputFiles/diff.K27.3d.USR.midfrag.bed") %>%
  filter(!(Chrom %in% c("Pt", "Mt"))) %>%
  filter(padj < 0.001) %>%
  filter(log2FC < -0.322) %>%
  dplyr::select(Chrom, Start, End, Length, Event, log2FC, pval, padj)

GRUnFilt <- makeGRangesFromDataFrame(UnFilt, keep.extra.columns = TRUE)

# Getting the gene lists --------------------------------------------------

GeneList <- function(x){
  y <- subsetByOverlaps(Genes, x, minoverlap = 150, ignore.strand = TRUE)
  y <- as.data.frame(unlist(y$GENEID))
  colnames(y) <- "Gene"
  left_join(y, TAIR10genes, by = "Gene")
}

GeneFilt <- GeneList(GRFilt)

GeneUnFilt <- GeneList(GRUnFilt)


# Venn diagram ------------------------------------------------------------
Venny <- list(Filt = GeneFilt$Gene,
              UnFilt = GeneUnFilt$Gene)
ggvenn(Venny)


# Looking for the lost genes ----------------------------------------------
Lost <- anti_join(GeneUnFilt, GeneFilt, by = "Gene")
