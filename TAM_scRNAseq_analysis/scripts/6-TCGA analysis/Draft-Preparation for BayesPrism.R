# ============================================================================
# Name: Preparation for BayesPrism
# Author: elopez
# Date: May-21-2024
# Description: 
# TODO: 
# ============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
# ============================================================================
# Section 1: DETERMINING MOST POLARIZED STATES
# ============================================================================
WuTAM_Pheno <-readRDS(file="D:/R_ScriptsPaper/Def_Objects/WuTAM_Pheno.rds")

# Extracting scores
scores <-colnames(WuTAM_Pheno@meta.data)[grepl("_score",colnames(WuTAM_Pheno@meta.data))]
# Polarized states
Polarized <- stringr::str_remove_all(scores,"_score")

Poles <- subset(WuTAM_Pheno, idents = Polarized)

# Convert the list to a named vector
Cells_ID <- names(Idents(Poles))
NewIdents <- paste0("TAM_",as.character(Idents(Poles)))

#Change idents in WuBRCA
WuBRCA_Refined<-readRDS(file = "D:/R_ScriptsPaper/Objects/WuBRCA_Refined.rds")
#First round of ident rename
WuBRCA_Refined <- SetIdent(WuBRCA_Refined, cells = Cells_ID, value = NewIdents)

# Remove undesired populations
WuBRCA_Refined <- subset(WuBRCA_Refined,
                         # Unlabeled TAMs, Doublets, and cycling cell types
                         idents = c("TAMs","Doublets","Cycling T CD8","Cycling TAMs"), invert=T)

WuBRCA_Refined <- RenameIdents(WuBRCA_Refined, 'Epithelial cells' = 'Tumor_cells')

#Save for BayesPrism
#saveRDS(WuBRCA_Refined, file="D:/BP_TME.rds")
