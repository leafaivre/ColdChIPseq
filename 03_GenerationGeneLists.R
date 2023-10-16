#This code generates the gene lists
#Uses the tracks generated in 01_GeneratingFilteredLists
#Author: Lea Faivre
#Date: 231012

# Libraries ---------------------------------------------------------------

library(plyr)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Athaliana.BioMart.plantsmart28)
library(readxl)
library(ChIPseeker)
library(xlsx)

Araport <- TxDb.Athaliana.BioMart.plantsmart28 

# Data import -------------------------------------------------------------
extraCols_FC <- c(Length = "numeric", Strand = "character", Event = "character", log2FC = "numeric", padj = "numeric")

DMR <- list(
  "K4_3h_Gain" = rtracklayer::import("Data/Test/K4hU.bed", format = "BED", 
                                  extraCols = extraCols_FC),
  "K4_3h_Loss" = rtracklayer::import("Data/Test/K4hD.bed", format = "BED", 
                                  extraCols = extraCols_FC),
  "K4_3d_Gain" = rtracklayer::import("Data/Test/K4dU.bed", format = "BED",
                                  extraCols = extraCols_FC),
  "K4_3d_Loss" = rtracklayer::import("Data/Test/K4dD.bed", format = "BED", 
                                  extraCols = extraCols_FC),
 "K27_3h_Gain" = rtracklayer::import("Data/Test/K27hU.bed", format = "BED", 
                                  extraCols = extraCols_FC),
  "K27_3h_Loss" = rtracklayer::import("Data/Test/K27hD.bed", format = "BED", 
                                  extraCols = extraCols_FC),
  "K27_3d_Gain" = rtracklayer::import("Data/Test/K27dU.bed", format = "BED", 
                                  extraCols = extraCols_FC),
  "K27_3d_Loss" = rtracklayer::import("Data/Test/K27dD.bed", format = "BED", 
                                  extraCols = extraCols_FC))

TAIR10genes <- read_excel("InputFiles/TAIR10genes.xlsx", 
                          sheet = "Name", col_types = c("text", "text"))

# Annotating the peaks ----------------------------------------------------
Anno <- lapply(DMR, annotatePeak, TxDb = Araport, level = "gene",
                 genomicAnnotationPriority = c("5UTR", "3UTR", "Exon", "Intron", 
                                               "Promoter", "Downstream", "Intergenic"),
                 tssRegion = c(-500, 500), verbose = FALSE)


AnnoDF <- list()
for (i in names(Anno)) {
  AnnoDF[[i]] <- data.frame(Anno[[i]]@anno)
}

# Keeping only the regions within genes -----------------------------------
AnnoDF <- lapply(AnnoDF, function(x){
  x %>% 
    filter(annotation != "Distal Intergenic",
           annotation != "Downstream (<=300bp)",
           annotation != "Downstream (<1kb)",
           annotation != "Downstream (1-2kb)",
           annotation != "Downstream (2-3kb)") %>%
    dplyr::select(seqnames, start, end, width, Event, log2FC, padj, annotation, 
                  geneId) 
})


# Correcting the missassignments ------------------------------------------
AnnoDFCorr <- lapply(AnnoDF, function(x){
  x %>%
    tidyr::separate(annotation, into = c("Anno", "gene"), sep = "\\(") %>%
    dplyr::mutate(gene = stringr::str_extract(gene, pattern = ".*(?=\\.)"),
                  Gene = ifelse(is.na(gene), geneId, gene)) %>%
    dplyr::select(-geneId, -gene) %>%
    left_join(TAIR10genes, by = "Gene")
})

DMGenes <- lapply(AnnoDFCorr, function(x){
  y <- x %>%
    dplyr::select(Gene)
  y <- unique(y)
})

##Number of uniques DM genes
sapply(DMGenes, nrow)

# Generating the files ----------------------------------------------------
for (i in names(AnnoDFCorr)) {
  if (nrow(AnnoDFCorr[[i]]) > 0) {
    write.xlsx(x = AnnoDFCorr[[i]],
               file = "Data/Test/DMGenes_Full.xlsx",
               sheetName = paste(i, sep = ""),
               row.names = FALSE,
               append = T)
  }
}
