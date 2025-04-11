# =============================================================================
# Name: TME ENRICHMENT
# Author: elopez
# Date: 19-08-2024
# Description: 
# TODO: 
# ============================================================================
#Intersecting signatures + rescued modules
BestGenes<-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_BestGenes.rds")
#
library(Seurat)
#Azizi dataset
#Infilt<-readRDS(file = "D:/R_ScriptsPaper/Objects/AzziBRCA_Refined.rds")
#Wu dataset
Infilt<-readRDS(file = "D:/R_ScriptsPaper/Objects/WuBRCA_Refined.rds")

Infilt<-subset(Infilt, idents = "Doublets", invert=T)
# ============================================================================
# Section 1: GSVA ENRICHMENT
# ============================================================================
datExpr<-GetAssayData(Infilt, assay = "RNA", slot = "data")
#Removing mitochondrial, Ribosomal and ORF genes
bad_patterns <- 'RP[0-9]|^MT|[0-9][0-9][0-9][0-9]|^RPL|^RPS|[0-9]orf[0-9]|^ATP'
a <- rownames(datExpr)[!grepl(bad_patterns, rownames(datExpr))]
datExpr<-datExpr[a,]
#Remove genes expressed in less than 5% of TAMs
datExpr<-datExpr[rowSums(datExpr == 0)/ncol(datExpr)<0.95,]
#Asses enrichment score with GSVA
library(GSVA)
Infilt_enrich<-as.data.frame(gsva(datExpr,BestGenes,method="gsva",parallel.sz=1))
#Agregando identidades a metadata
library(dplyr)
#A?adiendo scores enrichment a metadata
Tres<-as.data.frame(t(Infilt_enrich))
A<-colnames(Tres)
A<-paste(A, "score", sep = "_")
colnames(Tres)<-A
rm(datExpr)
#Preparando objeto seurat
#A?adiendo
Infilt_scores<-AddMetaData(Infilt,metada=Tres)
#saveRDS(Infilt_scores, file="D:/R_ScriptsPaper/Objects/Temporal/Infilt_scores_HSP.rds")
#saveRDS(Infilt_scores, file="D:/R_ScriptsPaper/Objects/Temporal/Azizi_Infilt_scores_HSP.rds")
#