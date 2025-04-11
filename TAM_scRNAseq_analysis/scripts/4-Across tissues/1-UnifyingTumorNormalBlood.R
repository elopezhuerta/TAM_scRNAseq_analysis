# =============================================================================
# Name: Unifying seurat objects
# Author: elopez
# Date: May-12-2024
# Description: 
# TODO: 
# ============================================================================
library(Seurat)
# ============================================================================
# Section 1: Blood from healthy volunteers (S.ObjA)
# ============================================================================
#Importando HBlood (Sangre de individuos sanos). Estan anotados
HBlood<-readRDS(file = "D:/R_ScriptsPaper/Objects/Wang_HBlood.rds")
#Seleccionando celulas anotadas como monocitos
#Quedandonos con monocytes
HBCells<-WhichCells(HBlood, idents = "Monocytes")
#Importando RawHBlood (sin normalizar, ni clusterizar ni anotar)
RawHBlood<-readRDS(file = "D:/R_Scripts/EssentialObjects/RawHBlood.rds")
#Solo las celulas que son monocitos
S.ObjA<-RawHBlood[,HBCells]
#Anadir columna de tejido, celltypes y tipo de paciente/individuo
S.ObjA@meta.data$Tissue<-"Blood"
S.ObjA@meta.data$TissueType<-"Blood: Healthy"
#Setting celltype and sample
S.ObjA@meta.data$celltypes<-"Monocytes"
S.ObjA@meta.data$lote<-S.ObjA$orig.ident
rm(HBlood,RawHBlood,HBCells)
#
# ============================================================================
# Section 2: Blood from BRCA patients (S.ObjB)
# ============================================================================
#Importando BrCBlood (Sangre de otra corte de pacientes con BrC). Estan anotados
BrCBlood<-readRDS(file = "D:/R_ScriptsPaper/Objects/Brechbuhl_BrCBlood.rds")
#Seleccionando celulas anotadas como monocitos
#Quedandonos con monocytes
BrCCells<-WhichCells(BrCBlood, idents = "Monocytes")
#Importando RawHBlood (sin normalizar, sin clusterisar o anotar)
RawBrCBlood<-readRDS(file = "D:/R_Scripts/EssentialObjects/RawBrCBlood.rds")
#Solo las celulas que son monocitos
S.ObjB<-RawBrCBlood[,BrCCells]
#Anadir columna de tejido, celltypes y tipo de paciente/individuo
S.ObjB@meta.data$Tissue<-"Blood"
S.ObjB@meta.data$TissueType<-"Blood: BRCA"
S.ObjB@meta.data$celltypes<-"Monocytes"
#Uniendo Batch 1 y 11 con 3 porque son poco abundantes
S.ObjB@meta.data$lote<-ifelse(S.ObjB$orig.ident%in%c("Batch 1","Batch 11",
                                                     "Batch 3","Batch 8,9,10"),
                              "Batch 1,3,8-11",
                              S.ObjB$orig.ident)
rm(BrCBlood,RawBrCBlood,BrCCells)
#
# ============================================================================
# Section 3: Non tumoral tissue from healthy volunteers (S.ObjC)
# ============================================================================
#Importando HNormal (Tejido de individuos sanos, operaciones est?ticas)
HNormal<-readRDS(file = "D:/R_ScriptsPaper/Objects/Twigger_HNormal.rds")
#Seleccionando celulas anotadas como monocitos macrofagos y DC (No hay DC)
HNCells<-WhichCells(HNormal, idents = "Tissue macrophages")
#Imporatndo RawHNormal (sin normalizar, clusterisar o anotar)
RawHNormal<-readRDS(file = "D:/R_Scripts/EssentialObjects/RawHNormal.rds")
#Solo las celulas que son monocitos
S.ObjC<-RawHNormal[,HNCells]
#Anadir columna de tejido, celltypes y tipo de paciente/individuo
#Ya tiene incluido el lote
S.ObjC$Tissue<-"Non-tumor tissue"
S.ObjC$TissueType<-"NT-tissue: Healthy"
S.ObjC$celltypes<-"Tissue macrophages"
#TIENE INCLUIDO EL LOTE
rm(HNormal,RawHNormal,HNCells)
# ============================================================================
# Section 4: Tumor, NT-tissue, blood from BRCA patients (TNB_MPC) (Azzi)
# ============================================================================
#Obtain only monocytes, TAMs and Tissue Mac
BRCA_Refined<-readRDS(file = "D:/R_ScriptsPaper/Objects/AzziBRCA_Refined.rds")
MPC<-WhichCells(BRCA_Refined, idents=c("TAMs","Tissue macrophages","Blood monocytes"))
#Importando tejido tumoral (Sangre, normal y tumoral)
AzziRaw<-readRDS(file = "D:/R_Scripts/EssentialObjects/AzziRaw.rds")
proj_name<-vector()
for (i in 1:length(AzziRaw)) {
  proj_name[i]<-unique(as.character(AzziRaw[[i]]$orig.ident))
}
#Merging
halb_seurat<-AzziRaw[c(-1)]
#Se anade peque?o identificador al ID de cada celula para que no este repetido (mismo batch)
merg_seurat<-merge(x=AzziRaw[[1]], y=halb_seurat, add.cell.ids = c(proj_name[1], proj_name[c(-1)]))

#Subseting MPC
TNB_MPC<-merg_seurat[,MPC]
#Save space
rm(BRCA_Refined, AzziRaw, merg_seurat,MPC,proj_name,halb_seurat)

#Adding tissue column
TNB_MPC$Tissue<-ifelse(grepl("BLOOD",TNB_MPC$orig.ident),
                       "Blood",
                       ifelse(grepl("NORMAL",TNB_MPC$orig.ident),
                              "Non-tumor tissue",
                              "Tumor"))
#Ya tiene datos de Tissue y celltypes. Falta agregar TissueType y lote
TNB_MPC$TissueType<-ifelse(TNB_MPC$Tissue=="Blood",
                           "Blood: BRCA",ifelse(TNB_MPC$Tissue=="Non-tumor tissue",
                                                "NT-tissue: BRCA",TNB_MPC$Tissue))

#Lote es igual que TissueType
TNB_MPC$lote<-TNB_MPC$TissueType

#Cambiando celltypes a (Monocytes, Tissue macrophages y TAM).
#Para monocytes falta especificar si son CD16 o CD14.
TNB_MPC$celltypes<-ifelse(TNB_MPC$Tissue=="Blood",
                          "Monocytes",ifelse(TNB_MPC$Tissue=="Non-tumor tissue",
                                               "Tissue macrophages","TAMs"))
#
# ============================================================================
# Section 5: Tumor TAMs (Wu_TAMs)
# ============================================================================
#Selected TAMs from Wu
TAMs<-readRDS(file="D:/R_ScriptsPaper/Objects/WuTAMs.rds")
#Raw Wu object
Wu_Raw<-readRDS(file = "D:/Munster/Muenster Workstation/huertaR/Objects/Seurat_V4/Wu_Seurat_V4.rds")
#Solo las celulas que son monocitos
Wu_TAMs<-Wu_Raw[,colnames(TAMs)]
Wu_TAMs$celltype_subset<-NULL
Wu_TAMs$celltype_major<-NULL
Wu_TAMs$celltype_minor<-NULL
Wu_TAMs$celltype_minor<-NULL
Wu_TAMs$celltype_minor<-NULL

Wu_TAMs$Tissue<-"Tumor"
Wu_TAMs$TissueType<-"Tumor"
Wu_TAMs$celltypes<-ifelse(colnames(Wu_TAMs)%in%WhichCells(TAMs, idents = "TAMs"),
                          "TAMs",
                          "Cycling TAMs")
#Lote es igual que sample
Wu_TAMs$lote<-TAMs$TAMsample
#Saving space
rm(TAMs,Wu_Raw)
# ============================================================================
# Section 6: Unifiying all the seurat objects
# ============================================================================
#GUARDANDO OBJETO CON monocitos y macrofagos de tejido no tumoral y sangre (de pacientes con BRCA y sanos)
RAW_MONOS_MAC<-merge(S.ObjA,list(S.ObjB, S.ObjC,TNB_MPC,Wu_TAMs))
#saveRDS(RAW_MONOS_MAC, file = "D:/R_ScriptsPaper/Objects/RAW_MONOS_MAC_Wu.rds")

Merg_list<-SplitObject(RAW_MONOS_MAC, split.by = "lote")

#No requieren QC. Pasar directo a normalizacion y batch correction
#Normalizacion
for (i in 1:length(Merg_list)) {
  Merg_list[[i]] <- NormalizeData(Merg_list[[i]], verbose = FALSE)
  Merg_list[[i]] <- FindVariableFeatures(Merg_list[[i]], selection.method = "vst", 
                                         nfeatures = 2000, verbose = FALSE)
}
#k.filter default cuando todos los archivos tienen mas de 200 celulas.
Merg.anchors <- FindIntegrationAnchors(object.list = Merg_list, dims = 1:30,
                                       k.filter = 150)
#IntegrateData returns new assay with 'batch-corrected' expression matrix for all cells, enabling them to be jointly analyzed.
Merg.anchors <- IntegrateData(anchorset = Merg.anchors, dims = 1:30)
# switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(Merg.anchors) <- "integrated"
Merg.anchors <- ScaleData(Merg.anchors, verbose = FALSE)
Merg.anchors <- RunPCA(Merg.anchors, npcs = 30, verbose = FALSE)
Merg.anchors <- RunUMAP(Merg.anchors, reduction = "pca", dims = 1:30)
Merg.anchors <- FindNeighbors(Merg.anchors, dims = 1:10)
#Higher resolition for more number of clusters
Merg.anchors <- FindClusters(Merg.anchors, resolution = 1.5)
#Fijando identidad
Idents(Merg.anchors)<-Merg.anchors$celltypes
# ============================================================================
# Section 7: Defying CD14 and CD16 monocytes
# ============================================================================
#Seleccionando monocitos
Monos <- subset(Merg.anchors,subset = celltypes=="Monocytes")
Monos <- RunTSNE(Monos, reduction = "pca", dims = 1:5)
Monos <- FindNeighbors(Monos, dims = 1:5)
Monos <- FindClusters(Monos, resolution = 0.2)

#Visualizando resultados
A<-DimPlot(Monos, reduction = "tsne", group.by = "seurat_clusters")
paleta<-RColorBrewer::brewer.pal(9, "YlOrRd")
B<-FeaturePlot(Monos, features = c("CD14","FCGR3A"), cols = paleta, reduction = "tsne")
A/B
#Cluster 2 son CD16
CD16<-WhichCells(Monos, expression=FCGR3A>=1)
Idents(Monos)<-ifelse(colnames(Monos)%in%CD16,"CD16 monocytes","CD14 monocytes")
#Cambiando identidades de AcrossTissues
Cells_ID<-names(Idents(Monos))
NewIdents<-as.character(Idents(Monos))
AcrossTissues <- SetIdent(Merg.anchors, cells = Cells_ID, value = NewIdents)

AcrossTissues$celltypes<-Idents(AcrossTissues)
#saveRDS(AcrossTissues, file = "D:/R_ScriptsPaper/Objects/AcrossTissues_Wu.rds")
#