# =============================================================================
# Name: WGCNA
# Author: huerta
# Date: May-03-2024
# Description:
# TODO: 
# =============================================================================
library(Seurat)
#Load seurat object with defined TAMs
#AZZI DATASET
S.Obj<-readRDS(file="D:/R_ScriptsPaper/Def_Objects/Azizi_TAMs.rds")

#WU DATASET
S.Obj<-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAMs.rds")

#Normalized counts. Genes as columns
Norm.Counts<-GetAssayData(S.Obj, assay = "RNA", slot = "data")
#This library is necessary
library(HelpersMG)
datExpr<-t(Norm.Counts)
#Structural genes are already removed
bad_patterns <- 'RP[0-9]|^MT|[0-9][0-9][0-9][0-9]|^RPL|^RPS|[0-9]orf[0-9]|^ATP'
a <- colnames(datExpr)[!grepl(bad_patterns, colnames(datExpr))]
datExpr<-datExpr[,a]
#Quitar genes que casi no estan expresados
datExpr<-datExpr[,colSums(datExpr == 0)/nrow(datExpr)<0.95]
#Wu: 8581 cells/8005 genes; Azizi:  2800 cells/3090 genes; softPower<-14

#Do not omit
options(stringsAsFactors = FALSE)
library(WGCNA)
# =============================================================================
# Section 1: Soft threshold (power)
# =============================================================================
#Signied matrix for positive correlation (not absolute values) and bicor for robust analysis in spite of outliers
#Picking softthreshold from the following range
powers<-c(c(1:10), seq(from = 12, to=30, by=2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", corFnc="bicor")
# Plot the results:
par(mfrow = c(1,2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", 
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], 
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")

#WU DATASET
#saveRDS(sft, file = "D:/R_ScriptsPaper/Def_Objects/Wu_SFT.rds")
#AZIZI DATASET
#saveRDS(sft, file = "D:/R_ScriptsPaper/Def_Objects/Azizi_SFT.rds")
# =============================================================================
# Section 2: Perform WGCNA
# =============================================================================
#Select sft
softPower<-12
#Module detection
#calculate the adjacency matrix
adj<-adjacency(datExpr, power = softPower, type = "signed", corFnc = "bicor")
#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
dissTOM<-1-TOMsimilarityFromExpr(datExpr, power = softPower, TOMType = "signed",  corType = "bicor")

#Module detection
#hierarchical clustering of the genes based on the TOM dissimilarity measure
library(flashClust)
geneTree<-flashClust(as.dist(dissTOM),method="average")

# Set the minimum module size
minModuleSize<-30
# Module ID using dynamic tree cut with recommended parameters (Intermediate deepsplit)
dynamicMods <- cutreeDynamic(dendro = geneTree, 
                             distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
#Change number for colors; note: The grey color is reserved for unassigned genes
dynamicColors<-labels2colors(dynamicMods)
#

#Merging similar modules
#Merge if modules are at least 80 similar 
MEDissThres<-0.20

# Call an automatic merging function
merge <- mergeCloseModules(as.matrix(datExpr), dynamicColors, cutHeight = MEDissThres)
# The merged module colors
mergedColors <- merge$colors
#To see what the merging did to our module colors
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Rename to dynamicColors
dynamicColors <- mergedColors

#Module similarity
#Calculates module eigengenes (1st principal component) of modules in a given single dataset.
MEList<- moduleEigengenes(as.matrix(datExpr), colors = dynamicColors)
datME<- MEList$eigengenes
plotEigengeneNetworks(datME, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2))
#

#Calculating gene memebership
datKME<-signedKME(datExpr, datME, outputColumnName="")
#Removing membership of grey genes
datKME<-datKME[,colnames(datKME)!="grey"]

#Extract every gene module
module_colors<-unique(dynamicColors)
if (length(module_colors)>=2){
  my_funct<-function(i){
    color<-module_colors[i]
    module<-colnames(datExpr)[which(dynamicColors==color)]
  }
  ModuleGenes<-lapply(c(1:length(module_colors)), FUN=my_funct)
  names(ModuleGenes)<-module_colors
} else {
  module<-colnames(datExpr)[which(dynamicColors==module_colors)]
  print(paste("Only one,",module_colors,"module",sep = " "))
}
str(ModuleGenes)

#WU DATASET
#saveRDS(datKME, file = "D:/R_ScriptsPaper/Def_Objects/Wu_datKME.rds")
#saveRDS(ModuleGenes, file = "D:/R_ScriptsPaper/Def_Objects/Wu_ModuleGenes.rds")

#AZIZI DATASET
#saveRDS(datKME, file = "D:/R_ScriptsPaper/Def_Objects/Azizi_datKME.rds")
#saveRDS(ModuleGenes, file = "D:/R_ScriptsPaper/Def_Objects/Azizi_ModuleGenes.rds")