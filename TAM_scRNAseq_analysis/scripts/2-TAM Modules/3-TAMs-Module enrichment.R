# =============================================================================
# Name: Module enrichment
# Author: Eric Lopez Huerta
# Date: May-08-2024
# Description:
# TODO: 
# =============================================================================
library(Seurat)
#CHOOSE DATASET TO WORK WITH
#Wu DATASET
TAMs<-readRDS(file="D:/R_ScriptsPaper/Def_Objects/WuTAMs.rds")

#AZZI DATASET
TAMs<-readRDS(file="D:/R_ScriptsPaper/Def_Objects/Azizi_TAMs.rds")

#Intersecting signatures + rescued modules
BestGenes<-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_BestGenes.rds")
#
# =============================================================================
# Section 1: Enrichment with BestGenes/Corrected Signatures
# =============================================================================
datExpr<-GetAssayData(TAMs, assay = "RNA", slot = "data")
#Removing mitochondrial, Ribosomal and ORF genes
bad_patterns <- 'RP[0-9]|^MT|[0-9][0-9][0-9][0-9]|^RPL|^RPS|[0-9]orf[0-9]|^ATP'
a <- rownames(datExpr)[!grepl(bad_patterns, rownames(datExpr))]
datExpr<-datExpr[a,]
#Remove genes expressed in less than 5% of TAMs
datExpr<-datExpr[rowSums(datExpr == 0)/ncol(datExpr)<0.95,]
#Asses enrichment score with GSVA
library(GSVA)
TAM_enrich<-as.data.frame(gsva(datExpr,BestGenes,method="gsva",parallel.sz=1))

#Adding enrichment data to seurat object
library(dplyr)
ES_DF<-as.data.frame(t(TAM_enrich))
A<-colnames(ES_DF)
scores<-paste(A, "score", sep = "_")
colnames(ES_DF)<-scores

#Choose location

#WU DATASET
#saveRDS(ES_DF, file = "D:/R_ScriptsPaper/Def_Objects/Wu_DF_ES_TAM.rds")

#DF with all the signatures in Azizi
#saveRDS(ES_DF, file = "D:/R_ScriptsPaper/Def_Objects/Azizi_DF_ES_TAM.rds")