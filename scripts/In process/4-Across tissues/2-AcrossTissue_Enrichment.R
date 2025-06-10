# =============================================================================
# Name: AcrossTissue_Enrichment
# Author: Eric Lopez Huerta
# Date: May-14-2024
# Description:
# TODO: 
# =============================================================================
library(Seurat)
#Intersecting signatures + rescued modules
BestGenes <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_BestGenes.rds")
# Import object with integrated monocytes, tissue macrophage and TAMs 
AcrossTissues<-readRDS(file="D:/R_ScriptsPaper/Objects/AcrossTissues_Wu.rds")
# =============================================================================
# Enrichment with corrected signatures
# =============================================================================
datExpr<-GetAssayData(AcrossTissues, assay = "RNA", slot = "data")
#Removing mitochondrial, Ribosomal and ORF genes
bad_patterns <- 'RP[0-9]|^MT|[0-9][0-9][0-9][0-9]|^RPL|^RPS|[0-9]orf[0-9]|^ATP'
a <- rownames(datExpr)[!grepl(bad_patterns, rownames(datExpr))]
datExpr<-datExpr[a,]
#Remove genes expressed in less than 5% of AcrossTissues
datExpr<-datExpr[rowSums(datExpr == 0)/ncol(datExpr)<0.95,]
#Asses enrichment score with GSVA
library(GSVA)
MonMacTAM_enrich<-as.data.frame(gsva(datExpr,BestGenes,method="gsva",parallel.sz=1))

#Adding enrichment data to seurat object
library(dplyr)
ES_DF<-as.data.frame(t(MonMacTAM_enrich))
A<-colnames(ES_DF)
scores<-paste(A, "score", sep = "_")
colnames(ES_DF)<-scores

#saveRDS(ES_DF, file = "D:/R_ScriptsPaper/Def_Objects/Tissues_ES_DF.rds")

# ============================================================================
# SHAPIRO TEST (NORMALITY TEST)
# ============================================================================
ES_DF <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Tissues_ES_DF.rds")
scores <- colnames(ES_DF)

S.Obj <-AddMetaData(AcrossTissues , metadata = ES_DF)

# Test by score per cell type
S.Obj_scores <-dplyr::group_split(S.Obj@meta.data[,c("TissueType",scores)], TissueType)

# Create a df indicating normality for each cell type per score
for (ii in seq_along(S.Obj_scores)) {
  print(ii)
  DFL <- S.Obj_scores[[ii]]
  #If cell type is > 5k, make sample because 5k is limit for shapiro test
  if (nrow(DFL)>5000) {
    DFL <- DFL[sample(nrow(DFL), 5000), ]
  }
  # Test for each module signature (per column), while removing TissueType column [,-1]
  shapiro_results <- apply(DFL[,-1], 2, shapiro.test)
  # Pval for each score
  shap_res <- vector()
  for (jj in 1:ncol(DFL[,-1])) {
    # If pval<0.05, data is not normally distributed
    shap_res[jj] <- shapiro_results[[jj]]$p.value > 0.05
  }
  result_df <- data.frame(names(shapiro_results),shap_res)
  colnames(result_df)<-c(as.character(unique(DFL$TissueType)),"Normal_distribution")
  print(result_df)
}
#
# ============================================================================
# KRUSKAL-WALLIS
# ============================================================================
library(dplyr)
library(ggpubr)
library(FSA)  # For Dunn's Test
# Were to store Dunn results
ns.comparisons <- list()
df.dunn.list <- list()

for (ii in seq_along(scores)) {
  # Create a df with only the desired columns.
  # Comparisons between different TissueTypes
  df <- S.Obj@meta.data[,c("TissueType",scores[ii])]
  # Homogenize names
  colnames(df) <- c("TissueType","ES")
  # Perform the Kruskal-Wallis Test
  kruskal_test <- kruskal.test(ES ~ TissueType, data = df)
  
  # If the Kruskal-Wallis test is significant, perform Dunn's Test for post-hoc analysis
  if (kruskal_test$p.value < 0.05) {
    dunn_test <- dunnTest(ES ~ TissueType, data = df, method = "bonferroni")
    dunn_res <- dunn_test$res
    # Save complete results in list of df
    df.dunn.list[[ii]] <- dunn_res
    # Rows were single.modules[ii] appear only once
    ns.pair <- dunn_res$Comparison[dunn_res$P.adj > 0.01]
    if (length(ns.pair)==0) {
      ns.pair <- "All significant"
    }
    ns.comparisons[[ii]] <- ns.pair
  } else {
    print(paste0("Kruskal-Wallis test was n.s.; no post-hoc test is performed for ",scores[ii]))
  }
}
# Summary of n.s. comparison with the respective population (single.modules[ii])
names(ns.comparisons)<-scores

ns.comparisons

# Whole results table
names(df.dunn.list) <- scores

# saveRDS(df.dunn.list, file = "D:/R_ScriptsPaper/Def_Objects/complete_tissue_comparisons.rds")
#