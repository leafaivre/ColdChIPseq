#This code generates the DEG lists from the featureCounts output.
#Author: Lea Faivre
#Date: 231024


# Libraries -----------------------------------------------------------------------------------
library(tidyverse)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(DEGreport)
library(readxl)
library(xlsx)

# Data import ---------------------------------------------------------------------------------

data <- read.table("InputFiles/withTE.txt", header = T, row.names = 1)

data <- data %>%
  select(-Chr, -Start, -End, -Strand, -Length)

colnames(data) <- c("N_1", "3h_1", "3d_1", 
                    "N_2", "3h_2", "3d_2", "N_3", "3h_3", "3d_3")

meta <- read.table("InputFiles/metadata_RNAseq.txt", header = T, row.names = 1)
meta <- meta %>% 
  mutate(Rep = as.factor(Rep))


# Creating the DESeq2 object ------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData = data, 
                              colData = meta, 
                              design = ~ Rep + Condition)


# Normalization -------------------------------------------------------------------------------

dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = TRUE)

# PCA -----------------------------------------------------------------------------------------

rld <- rlog(dds, blind = TRUE)

pcaData <- plotPCA(rld, intgroup = c("Condition", "Rep"), returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = Condition, shape = Rep)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("blue4", "cornflowerblue", "darkgreen")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed() +
  theme_bw()
#ggsave("Plots/PCA_RNAseq.tiff", units = "in", width = 4, height = 4, dpi = 300, compression = 'lzw')


# DE analysis ---------------------------------------------------------------------------------

# Preliminary checks ___________________________________________________________
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

# DE analysis _______________________________________________________________

dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ Rep + Condition)

dds <- DESeq(dds)

Contrasts <- list(
  Nvs3h = c("Condition", "3h", "N"),
  Nvs3d = c("Condition", "3d", "N"))

Results <- lapply(Contrasts, function(x){
  y <- results(dds, contrast = x, alpha = 0.05)
  y_shr <- lfcShrink(dds, contrast = x, res = y, type = "ashr")
})


# Extracting the significant DEG --------------------------------------------------------------

Res_all <- lapply(Results, function(x){
  x %>%
    data.frame() %>%
    rownames_to_column(var = "Gene") %>%
    as_tibble()})

Sig_Res <- lapply(Res_all, function(x){
  x %>% 
    filter(padj < 0.05 & abs(log2FoldChange) >= 1)
})


# Writing the files ---------------------------------------------------------------------------
write.xlsx(x = as.data.frame(Res_all[["Nvs3h"]]),
           file = "Data/RNAseq_Nvs3h.xlsx",
           row.names = FALSE)

write.xlsx(x = as.data.frame(Res_all[["Nvs3d"]]),
           file = "Data/RNAseq_Nvs3d.xlsx",
           row.names = FALSE)

for (i in names(Sig_Res)) {
  write.xlsx(x = as.data.frame(Sig_Res[[i]]),
             file = "Data/DEGenes.xlsx",
             sheetName = paste(i, sep = ""),
             row.names = FALSE,
             append = T)
}
