# ============================================================================
# Name: M1_M2 Enrichment 
# Author: elopez
# Date: Sep-27-2024
# Description: ES of curated M1 M2 signatures in TAMs and TME.
# TODO: 
# ============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
#Color palette
reds<-RColorBrewer::brewer.pal(9, "YlOrRd")
single.modules <- c("FCN1","ISG15","CXCL9","FOSB","HLA_II","HSPs","C1QA","APOE")
#Importing signature to test
file_path<-"D:/Antecedentes y teoría/Macrofagos/M1_M2/1-Tabla definitica Referencias M1 M2.xlsx"
pre_M1 <- readxl::read_excel(file_path, 
                              sheet = "M1 markers",
                              range = "A1:A22", col_names = T)

pre_M2 <- readxl::read_excel(file_path, 
                             sheet = "M2 markers",
                             range = "A1:C34", col_names = T)
# As list
M1_M2 <- list(
  M1_Curated = pre_M1$`M1 marker`,
  M2_Curated = pre_M2$`M2 marker`[!is.na(pre_M2$`Reference (DOI)`)]
)

# Importing object
WuBRCA<-readRDS(file = "D:/R_ScriptsPaper/Objects/WuBRCA_Refined.rds")
# Removing doublets for visualization
WuBRCA<-subset(WuBRCA, idents = "Doublets", invert=T)

# Importing TAMs from Wu dataset
TAM_Modules<-readRDS(file = "D:/R_ScriptsPaper/Objects/Temporal/WuTAM_PhenoHSP.rds")
#
# ============================================================================
# Enrichment in TME
# ============================================================================
datExpr<-GetAssayData(WuBRCA, assay = "RNA", slot = "data")
#Removing mitochondrial, Ribosomal and ORF genes
bad_patterns <- 'RP[0-9]|^MT|[0-9][0-9][0-9][0-9]|^RPL|^RPS|[0-9]orf[0-9]|^ATP'
a <- rownames(datExpr)[!grepl(bad_patterns, rownames(datExpr))]
datExpr<-datExpr[a,]
#Remove genes expressed in less than 5% of TAMs
datExpr<-datExpr[rowSums(datExpr == 0)/ncol(datExpr)<0.95,]
#Asses enrichment score with GSVA
library(GSVA)
M1M2.Enrichment<-as.data.frame(gsva(datExpr,M1_M2,method="gsva",parallel.sz=1))
# Turn into DF
library(dplyr)
DF_enrich<-as.data.frame(t(M1M2.Enrichment))
# Modifying colnames
colnames(DF_enrich)<-paste0(colnames(DF_enrich), "_score")
rm(datExpr)
# Saving results
saveRDS(DF_enrich, file = "D:/R_ScriptsPaper/Objects/Temporal/ES_M1M2_WuBRCA.rds")

#TEMPORAL
M1M2.WuBRCA<-AddModuleScore(WuBRCA, features = M1_M2,name = "M", assay = "RNA")
colnames(M1M2.WuBRCA@meta.data)<-ifelse(colnames(M1M2.WuBRCA@meta.data)%in%c("M1","M2"),
                                        paste(colnames(M1M2.WuBRCA@meta.data),"score", sep = " "),
                                        colnames(M1M2.WuBRCA@meta.data))
# Save 
saveRDS(M1M2.WuBRCA, file = "D:/R_ScriptsPaper/Objects/Temporal/ModuleScoreM1M2_WuBRCA.rds")
#
# ============================================================================
# Enrichment in TAMs
# ============================================================================

datExpr<-GetAssayData(TAM_Modules, assay = "RNA", slot = "data")
#Removing mitochondrial, Ribosomal and ORF genes
bad_patterns <- 'RP[0-9]|^MT|[0-9][0-9][0-9][0-9]|^RPL|^RPS|[0-9]orf[0-9]|^ATP'
a <- rownames(datExpr)[!grepl(bad_patterns, rownames(datExpr))]
datExpr<-datExpr[a,]
#Remove genes expressed in less than 5% of TAMs
datExpr<-datExpr[rowSums(datExpr == 0)/ncol(datExpr)<0.95,]
#Asses enrichment score with GSVA
library(GSVA)
M1M2.Enrichment<-as.data.frame(gsva(datExpr,M1_M2,method="gsva",parallel.sz=1))
# Turn into DF
library(dplyr)
DF_enrich<-as.data.frame(t(M1M2.Enrichment))
# Modifying colnames
colnames(DF_enrich)<-paste0(colnames(DF_enrich), "_score")
rm(datExpr)
# Saving results
saveRDS(DF_enrich, file = "D:/R_ScriptsPaper/Objects/Temporal/ES_M1M2_WuTAMs.rds")

# ============================================================================
# DIFF EXPRESSION M1 M2
# ============================================================================
#Import table of DEG
DEG_Polos <- readRDS(file = "D:/R_ScriptsPaper/Objects/Temporal/DEG_POLOS_HSP.rds")

# Define M1_M2 genes as a character vector
M1_M2_genes <- as.character(unlist(M1_M2))

# Filter and arrange the dataframe
DEGM1M2 <- DEG_Polos %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  filter(avg_log2FC > 0,                               # Only positive avg_log2FC values
         gene %in% M1_M2_genes,                        # Filter by genes in M1_M2
         p_val_adj < 0.05)                             # Filter by adjusted p-value < 0.05

NotSig_M1M2<-M1_M2_genes[!M1_M2_genes%in%unique(DEGM1M2$gene)]

# Display the result for genes with p_val_adj > 0.05
print(NotSig_M1M2)

#