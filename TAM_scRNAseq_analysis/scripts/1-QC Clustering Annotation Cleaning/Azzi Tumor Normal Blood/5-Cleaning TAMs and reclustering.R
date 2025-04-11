# =============================================================================
# Name: Cleaning and reclustering of TAMs
# Author: elopez
# Date: Apr-21-2024
# Description: 

# TODO: 
# =============================================================================
library(Seurat)
library(RColorBrewer)
library(ggplot2)
paleta<-brewer.pal(9, "YlOrRd")
Annotated<-readRDS("D:/R_ScriptsPaper/Objects/AzziAnnotated.rds")
# ============================================================================
# Section 1: 
# ============================================================================
#Label sparse cells as doublets interactively
Sub_Annot<-subset(Annotated, subset=Annotation %in% c("TAMs"))
#Select sparse cells interactively
#This are possible duplets (cells in lymphoid clusters)  
T_Doublets<-CellSelector(plot = DimPlot(Sub_Annot))
B_Doublets<-CellSelector(plot = DimPlot(Sub_Annot))
#
# ============================================================================
# Section 2: FIRST RECLUSTERING
# ============================================================================
#Select true TAMs
TAMs.ID<-WhichCells(Annotated, idents="TAMs")
#Remove duplets
TAMs.ID<-TAMs.ID[!TAMs.ID%in%c(T_Doublets,B_Doublets)]
#List of Seurat Objects with no normalization, clustering or annotation 
my_seurat<-readRDS(file = "D:/R_Scripts/EssentialObjects/AzziRaw.rds")

proj_name<-vector()
for (i in 1:length(my_seurat)) {
  proj_name[i]<-unique(as.character(my_seurat[[i]]$orig.ident))
}

#Merging
halb_seurat<-my_seurat[c(-1)]
#Se a?ade peque?o identificador al ID de cada celula para que no este repetido (mismo batch)
merg_seurat<-merge(x=my_seurat[[1]], y=halb_seurat, add.cell.ids = c(proj_name[1], proj_name[c(-1)]))

#Select only TAMs
S.Obj<-merg_seurat[,TAMs.ID]

#Adding sample column with summarized information
S.Obj$Sample<-stringr::str_sub(S.Obj$orig.ident,start = 0,end = 4)
#Merging BC02 into BC08 to augment cell size withing sample.
S.Obj$Sample<-ifelse(S.Obj$Sample%in%c("BC02","BC08"),"BC02_08",S.Obj$Sample)

Merg_list<-SplitObject(S.Obj, split.by = "Sample")
for (i in 1:length(Merg_list)) {
  Merg_list[[i]] <- NormalizeData(Merg_list[[i]], verbose = FALSE)
  Merg_list[[i]] <- FindVariableFeatures(Merg_list[[i]], selection.method = "vst", 
                                         nfeatures = 2000, verbose = FALSE)
}
#k.filter default
Merg.anchors <- FindIntegrationAnchors(object.list = Merg_list, dims = 1:30)
#IntegrateData returns new assay with 'batch-corrected' expression matrix for all cells, enabling them to be jointly analyzed.
Merg.anchors <- IntegrateData(anchorset = Merg.anchors, dims = 1:30)
# switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(Merg.anchors) <- "integrated"
Merg.anchors <- ScaleData(Merg.anchors, verbose = FALSE)
Merg.anchors <- RunPCA(Merg.anchors, npcs = 30, verbose = FALSE)
Merg.anchors <- RunUMAP(Merg.anchors, reduction = "pca", dims = 1:30)
Merg.anchors <- FindNeighbors(Merg.anchors, dims = 1:10)
Merg.anchors <- FindClusters(Merg.anchors, resolution = 0.5)
#
# ============================================================================
# Section 3: REFINMENT OF BRCA_Refined BASED ON CLUSTER AND RESPECTIVE MARKERS
# ============================================================================
#Marker finding for each cluster
ClustMarkers<-FindAllMarkers(Merg.anchors,only.pos = TRUE,min.pct = 0.25, min.cells.feature = 15)
library(dplyr)
#List of first markers for each cluster, considering 25 genes per cluster
DF_Mark<-ClustMarkers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)
#List with first 10 markers for each cluster
SIGNATURES<-list()
#Vector with names for each SIGNATURE and PATH
SigName<-vector()
#Every marker for each cluster
PATHS<-list()
for(i in 1:length(unique(DF_Mark$cluster))){
  PATHS[[i]]<-DF_Mark$gene[DF_Mark$cluster==c(i-1)]
  #First 10 markers for enrichment
  SIGNATURES[[i]]<-DF_Mark$gene[DF_Mark$cluster==c(i-1)]
  SIGNATURES[[i]]<-SIGNATURES[[i]][1:10]
  SigName[i]<-SIGNATURES[[i]][1]
}
names(SIGNATURES)<-SigName
names(PATHS)<-SigName
#
# IDENTIFIYING TAMS ENRICHED WITH SIGNATURES NOT BELOGING TO TAMS
#Using signature with 10 genes
SigKey<-paste(names(SIGNATURES),"score", sep = " ")
SigKey<-stringr::str_replace(SigKey,"-"," ")
#Enrichment with z-scoring. (Tirosh et al. 2017)
#single cell has few genes. Use lower number for ctrl parameter in AddModuleScore() funciton
Sub_Infilt_score<-AddModuleScore(Merg.anchors, features = SIGNATURES, ctrl = 50, name = SigKey, assay = "RNA")
#Modifying names
title<-stringr::str_replace_all(colnames(Sub_Infilt_score@meta.data), "e[0-9][0-9]|e[0-9]", "e")
colnames(Sub_Infilt_score@meta.data)<-title
scores<-colnames(Sub_Infilt_score@meta.data)[grep("score",colnames(Sub_Infilt_score@meta.data))]
#Cells expressing each signature in SIGNATURE
library(RColorBrewer)
library(ggplot2)
paleta<-brewer.pal(9, "YlOrRd")

#Visuals
DimPlot(Sub_Infilt_score, pt.size = 0.9, label = T)+NoLegend()+
  xlab(NULL)+ ylab(NULL)+
  ggtitle("Myeloid cells")

FeaturePlot(Sub_Infilt_score, features = scores, cols=paleta,pt.size = 0.8)+
  xlab(NULL)+ ylab(NULL)
#
# SIGNATURE ENRICHMENT IN THE MICROENVIRONMENT
Annotated_Score<-AddModuleScore(Annotated, features = SIGNATURES, ctrl = 19, name = SigKey)
#Modifying names
title<-stringr::str_replace_all(colnames(Annotated_Score@meta.data), "e[0-9][0-9]|e[0-9]", "e")
colnames(Annotated_Score@meta.data)<-title
scores<-colnames(Annotated_Score@meta.data)[grep("score",colnames(Annotated_Score@meta.data))]
#Cells expressing each signature in SIGNATURE
curvas<-function(i){
  iPlot<-RidgePlot(Annotated_Score, features = scores[i])+
    #El texto orizontal que etiqueta cada curva
    theme(axis.text = element_text(size = rel(0.6)),
          #El texto en vertical
          axis.title.y = element_blank(),
          #El titulo
          plot.title = element_text(size=rel(0.7),
                                    hjust = 0.5))+
    #Sin la leyenda
    NoLegend()
}
lista_Ridge<-lapply(c(1:length(scores)),FUN = "curvas")
gridExtra::grid.arrange(grobs = lista_Ridge, ncol = 4)
#
# PATHWAYS REPRESENTED BY EACH SIGNATURE PATH
library(org.Hs.eg.db)
library(limma)
library(dplyr)
library(tibble)
path_enrich<-function(i){
  color<-as.data.frame(PATHS[[i]])
  color$EntrezID <- mapIds(org.Hs.eg.db,keys=color[,1],column="ENTREZID",keytype="SYMBOL",multiVals="first")
  #enrichment with GO and KEGG databases
  goRest<-goana(color$EntrezID, species="Hs")
  k<-kegga(color$EntrezID, species="Hs")
  #Sorting pathways with pval<0.01 in both databases (GO/KEGG)
  goSort<-goRest%>%
    filter(Ont%in%"BP"&P.DE<0.01)%>%
    arrange(P.DE)
  KEGGSort<-k%>%
    filter(P.DE<0.01)%>%
    arrange(P.DE)
  KEGGSort$Ont<-rep("KEGG",nrow(KEGGSort))
  KEGGSort<-KEGGSort[,colnames(KEGGSort)[c(1,5,2:4)]]
  colnames(KEGGSort)<-colnames(goSort)
  pathways<-rbind(goSort,KEGGSort)
  rownames_to_column(pathways, var = "ID")
}
pathways<-lapply(c(1:length(PATHS)),FUN = "path_enrich")
#Choosing top 6 pathways of each database (GO/KEGG)
path_plots<-vector("list", length(pathways))
for (i in 1:length(pathways)){
  path<-pathways[[i]]
  #Removing pathways with more than 1000 genes beacause they are not descriptive enough
  path_plots[[i]]<-path%>%filter(N<900)%>%
    #Top 6 pathways of GO and KEGG
    group_by(Ont)%>%
    slice_min(n = 5, order_by = P.DE)
}
names(path_plots)<-names(PATHS)

#Concatenating every dataframe (DF) into a list of DF
list_DF<-function(i){
  #Creating vector with name of every Cell cluster (Module)
  neun_col<-rep(names(path_plots[i]),nrow(path_plots[[i]]))
  #Adding column with names of their respective cell cluster (Module)
  path_plots[[i]]$"Module"<-neun_col
  DFL<-as.data.frame(path_plots[[i]])
}
DF_list<-lapply(c(1:length(path_plots)),FUN = "list_DF")
#Concatenating every DF into single unified DF
D<-do.call("rbind",DF_list)

#Separating GO results from KEGG
library(dplyr)
DL<-group_split(D, Ont)
library(ggplot2)
#Choose database to creat graph
base<-c("Gene Ontology", "KEGG")
DF<-as.data.frame(DL[[1]])
#Setting desired order for pathway (Term) and cell cluster (Module)
Pathway<-factor(DF$Term, levels=unique(DF$Term))
Cell_cluster<-factor(DF$Module, levels=unique(DF$Module))
DF$Term<-Pathway
DF$Module<-Cell_cluster
#Calculating ratio of genes participating (DF$DE) in the complete list of genes (DF$N) in every pathway
DF$GeneRatio<-DF$DE/DF$N
#Colour palette
library(RColorBrewer)
paleta<-brewer.pal(9, "YlOrRd")
ggplot(DF, aes(x= Cell_cluster, y= Pathway, fill=GeneRatio))+
  geom_tile()+
  #setting white background
  theme_classic()+
  #Choosing colour and title
  scale_fill_gradientn(colours = paleta)+ggtitle(paste("Pathways enriched according to",base[1]))+
  #x thick marks orientation and vertical dashlines 
  theme(axis.text.x=element_text(angle = 30, hjust = 1, vjust = 1),
        panel.border = element_rect(fill=NA),
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank())

#
#DECISION
#REMOVE (IF NECCESSARY):
#CCL5 (lymphoid NK/T/B)
#Remove Clusters 9 and 11
#Labeling cells with kmeans, dividing in 4 groups (Low, Intermediate 1, Intermediate 2, High)
cell_contamination<-c("CCL5.score")
remove_cells<-list()
for (i in 1:length(cell_contamination)){
  #for each signature, separate cells in 4 groups with kmeans
  cont_sub<-Sub_Infilt_score@meta.data[,cell_contamination[i]]
  sep<-kmeans(cont_sub, 4)
  #find cluster with the higher value
  cl<-names(sep$centers[which.max(sep$centers), ])
  remove_cells[[i]]<-rownames(Sub_Infilt_score@meta.data[sep$cluster==cl,])
}
#adding doublets previously selecte
remove_cells[[2]]<-c(T_Doublets,B_Doublets)
#Renaming missing cells as doublets
remove_cells[[3]]<-TAMs.ID[!TAMs.ID%in%colnames(merg_seurat)]

#Removing cells
to_remove<-unique(as.character(unlist(remove_cells)))

#Renombrando como lymphoid myeloid
BRCA_Refined <- SetIdent(Annotated, cells = to_remove, value = "Doublets")
BRCA_Refined$celltypes<-Idents(BRCA_Refined)
BRCA_Refined$Annotation<-NULL

#Saving
saveRDS(BRCA_Refined, file = "D:/R_ScriptsPaper/Objects/AzziBRCA_Refined.rds")
#

# ============================================================================
# Section 4: RECLUSTERING MYELOID CELLS
# ============================================================================
BRCA_Refined<-readRDS(file = "D:/R_ScriptsPaper/Objects/AzziBRCA_Refined.rds")
Sub_Annot<-subset(BRCA_Refined,idents = c("TAMs","DC populations","Neutrophils","Mast cells","pDC"))
#Only cells from tumor
tumor_assoc<-Sub_Annot$orig.ident[grep("TUMOR",Sub_Annot$orig.ident)]
#Tumor infiltrating myeloid
Myeloid<-subset(Sub_Annot, subset=orig.ident %in% tumor_assoc)

#Clustering TAM, monocytes and DC, mast cells and pDC with tsne
Myeloid <- RunTSNE(Myeloid, reduction = "pca", dims = 1:30)
Myeloid <- FindNeighbors(Myeloid, dims = 1:30)
Myeloid <- FindClusters(Myeloid)

Idents(Myeloid)<-Myeloid$celltypes
saveRDS(Myeloid, file = "D:/R_ScriptsPaper/Objects/Myeloid.rds")
#
#Myeloid markers analysis
#Mac/Mo: 
feats<-c("CD68","CD163","CD14","FCGR3A", 
         #DC      
         "THBD","IRF8","ITGAX","FCER1A","CD1C","IRF4",
         #pDC     
         "LILRA4","GZMB", "CXCR3", 
         #Mastcells 
         "TPSAB1","CPA3", 
         #Neutrophils 
         "FCGR3B","S100A9")

#Displaying myeloid markers with custom order
plotear<-function(i){
  ip<-FeaturePlot(Myeloid, features = feats[i], cols=paleta, reduction = "tsne",pt.size = 0.8)+
    xlab(NULL)+ ylab(NULL)
}
lista_heat<-lapply(c(1:length(feats)),FUN = "plotear")
ggpubr::ggarrange(plotlist = lista_heat, ncol=7, nrow=3, common.legend = TRUE)
#
# ============================================================================
# Section 5: SECOND TAM RECLUSTERING
# ============================================================================
#List of Seurat Objects with no normalization, clustering or annotation 
my_seurat<-readRDS(file = "D:/R_Scripts/EssentialObjects/AzziRaw.rds")

proj_name<-vector()
for (i in 1:length(my_seurat)) {
  proj_name[i]<-unique(as.character(my_seurat[[i]]$orig.ident))
}
#Merging
halb_seurat<-my_seurat[c(-1)]
#Se anade pequeno identificador al ID de cada celula para que no este repetido (mismo batch)
merg_seurat<-merge(x=my_seurat[[1]], y=halb_seurat, add.cell.ids = c(proj_name[1], proj_name[c(-1)]))

#Select only TAMs from BRCA_Refined
BRCA_Refined<-readRDS(file = "D:/R_ScriptsPaper/Objects/AzziBRCA_Refined.rds")
TAMs.ID<-WhichCells(BRCA_Refined, idents = "TAMs")
S.Obj<-merg_seurat[,TAMs.ID]

#Adding sample column with summarized information
S.Obj$Sample<-stringr::str_sub(S.Obj$orig.ident,start = 0,end = 4)
#Merging BC02 into BC08 to augment cell size withing sample.
S.Obj$Sample<-ifelse(S.Obj$Sample%in%c("BC02","BC08"),"BC02_08",S.Obj$Sample)

Merg_list<-SplitObject(S.Obj, split.by = "Sample")
for (i in 1:length(Merg_list)) {
  Merg_list[[i]] <- NormalizeData(Merg_list[[i]], verbose = FALSE)
  Merg_list[[i]] <- FindVariableFeatures(Merg_list[[i]], selection.method = "vst", 
                                         nfeatures = 2000, verbose = FALSE)
}
#k.filter default
Merg.anchors <- FindIntegrationAnchors(object.list = Merg_list, dims = 1:30)
#IntegrateData returns new assay with 'batch-corrected' expression matrix for all cells, enabling them to be jointly analyzed.
Merg.anchors <- IntegrateData(anchorset = Merg.anchors, dims = 1:30)
# switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(Merg.anchors) <- "integrated"
Merg.anchors <- ScaleData(Merg.anchors, verbose = FALSE)
Merg.anchors <- RunPCA(Merg.anchors, npcs = 30, verbose = FALSE)
Merg.anchors <- RunUMAP(Merg.anchors, reduction = "pca", dims = 1:30)
Merg.anchors <- FindNeighbors(Merg.anchors, dims = 1:10)
Merg.anchors <- FindClusters(Merg.anchors, resolution = 0.5)
#Saving
saveRDS(Merg.anchors, file = "D:/R_ScriptsPaper/Def_Objects/Azizi_TAMs.rds")
