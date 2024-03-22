library(BiocManager)
library(DiffBind)
library(tidyverse)

setwd("/scratch/leafaivre/ColdChIP_new/results")

samples<- read.csv("~/ColdChIP_new/ChIPQC/K27_USR.csv")

ObjUS <- dba(sampleSheet = samples)

ObjUS <- dba.count(ObjUS, score = DBA_SCORE_READS, bUseSummarizeOverlaps = TRUE, minOverlap = 1)

K27 <- dba.contrast(ObjUS, categories = DBA_CONDITION, minMembers = 2)
save(K27, file = "DiffBind/K27")


samples4<- read.csv("~/ColdChIP_new/ChIPQC/K4_USR.csv")

Obj4 <- dba(sampleSheet = samples4)

Obj4 <- dba.count(Obj4, score = DBA_SCORE_READS, bUseSummarizeOverlaps = TRUE, minOverlap = 1)

K4 <- dba.contrast(Obj4, categories = DBA_CONDITION, minMembers = 2)
save(K4, file = "DiffBind/K4")





