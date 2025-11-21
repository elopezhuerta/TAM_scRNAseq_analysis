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
Annotated<-readRDS("D:/R_ScriptsPaper/Objects/WuBRCA_Annotated.rds")
# ============================================================================
# Section 1: Interactive selection
# ============================================================================
#Label sparse cells as doublets interactively
Sub_Annot<-subset(Annotated, subset=celltypes %in% c("TAMs", "Cycling TAMs"))
#Select TAMs clustered together
TAMs.ID<-CellSelector(plot = DimPlot(Sub_Annot))
#These are possible duplets (cells in lymphoid clusters)  
Doublets<-colnames(Sub_Annot)[!colnames(Sub_Annot)%in%TAMs.ID]
#
# ============================================================================
# Section 2: FIRST RECLUSTERING
# ============================================================================
#List of Seurat Objects with no normalization, clustering or annotation 
my_seurat<-readRDS(file = "D:/Munster/Muenster Workstation/huertaR/Objects/Seurat_V4/Wu_Seurat_V4.rds")
my_seurat@meta.data[,c("celltype_subset","celltype_minor","celltype_major")]<-NULL

#Select only TAMs
S.Obj<-my_seurat[,TAMs.ID]
rm(my_seurat)

#Some samples have less than 200 cells. We have to group them in pairs
samp.size<-sort(table(S.Obj$orig.ident),decreasing=T)
samp.name<-names(samp.size)[samp.size<=200]
#Half the samples
first_half<-samp.name[1:(length(samp.name)/2)]
second_half<-rev(samp.name)[1:(length(samp.name)/2)]
#New name in pairs
fir_sec<-paste(first_half,second_half,sep = "_")
#New sample grouping
S.Obj$TAMsample<-as.character(S.Obj$orig.ident)
for (i in 1:length(first_half)) {
  S.Obj$TAMsample<-ifelse(as.character(S.Obj$TAMsample)%in%c(first_half[i],second_half[i]),
                         fir_sec[i],
                         as.character(S.Obj$TAMsample))
}
#Batch correction
#1-Split by sample
Merg_list <- SplitObject(S.Obj, split.by = "TAMsample")
#2-Normalize cells from each individual sample
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
Merg.anchors <- FindClusters(Merg.anchors, resolution = 1)
#
# ============================================================================
# Section 3: REFINMENT OF WuBRCA_Annotated BASED ON CLUSTER AND RESPECTIVE MARKERS
# ============================================================================
#Marker finding for each cluster
ClustMarkers<-FindAllMarkers(Merg.anchors,only.pos = TRUE, assay = "RNA",
                             min.pct = 0.25, min.cells.feature = 15)
library(dplyr)
#List of first markers for each cluster, considering 25 genes per cluster
DF_Mark<-ClustMarkers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
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
Merg.anchors.scores<-AddModuleScore(Merg.anchors, features = SIGNATURES, ctrl = 20, name = SigKey, assay = "RNA")
#Modifying names
title<-stringr::str_replace_all(colnames(Merg.anchors.scores@meta.data), "e[0-9][0-9]|e[0-9]", "e")
colnames(Merg.anchors.scores@meta.data)<-title
scores<-colnames(Merg.anchors.scores@meta.data)[grep("score",colnames(Merg.anchors.scores@meta.data))]
#Cells expressing each signature in SIGNATURE
library(RColorBrewer)
library(ggplot2)
paleta<-brewer.pal(9, "YlOrRd")

#Visuals
DimPlot(Merg.anchors.scores, pt.size = 0.9, label = T)+NoLegend()+
  xlab(NULL)+ ylab(NULL)+
  ggtitle("Myeloid cells")

FeaturePlot(Merg.anchors.scores, features = scores, cols=paleta,pt.size = 0.8)+
  xlab(NULL)+ ylab(NULL)


# SIGNATURE ENRICHMENT IN THE MICROENVIRONMENT
Annotated_Score<-AddModuleScore(Annotated, features = SIGNATURES, 
                                ctrl = 19,name = SigKey,assay = "RNA")
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
# ============================================================================
# Section 4: PATHWAYS REPRESENTED BY EACH SIGNATURE PATH
# ============================================================================
library(org.Hs.eg.db)
library(limma)
library(dplyr)
library(tibble)

#Descargando lista de vias para KEGG
PathDF<-getKEGGPathwayNames(species="hsa", remove.qualifier = T)
colnames(PathDF)<-c("ID","Pathway")
#PathDF<-column_to_rownames(PathDF, var = "PathwayID")
path_enrich<-function(i){
  print(i)
  color<-as.data.frame(PATHS[[i]])
  color$EntrezID <- mapIds(org.Hs.eg.db,keys=color[,1],
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
  #enrichment with GO and KEGG databases
  goRest<-goana(color$EntrezID, species="Hs")
  k<-kegga(color$EntrezID, species="Hs")
  #A?adiendo nombre largo a KEGG
  k$ID<-stringr::str_remove_all(rownames(k),"path:")
  k$Pathway<-NULL
  k<-left_join(k,PathDF, by="ID")
  #Dejandola en rownames (donde estaban)
  k<-column_to_rownames(k, var="ID")
  #Sorting pathways with pval<0.01 in both databases (GO/KEGG)
  goSort<-goRest%>%
    filter(Ont%in%"BP"&P.DE<0.01)%>%
    arrange(P.DE)
  KEGGSort<-k%>%
    filter(P.DE<0.01)%>%
    arrange(P.DE)
  KEGGSort$Ont<-rep("KEGG",nrow(KEGGSort))
  #Re ordenando columnas de KEGGSort
  KEGGSort<-KEGGSort[,order(match(colnames(KEGGSort),c("Pathway","Ont","N","DE","P.DE")))]
  colnames(KEGGSort)<-colnames(goSort)
  pathways<-rbind(goSort,KEGGSort)
  rownames_to_column(pathways, var = "ID")
}
pathways<-lapply(c(1:length(PATHS)),FUN = "path_enrich")
#Choosing top 5 pathways of each database (GO/KEGG)
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
DL<-group_split(D, Ont)
library(ggplot2)
#Colour palette
paleta<-RColorBrewer::brewer.pal(9, "Greens")[2:9]
#Choose database to creat graph
base<-c("GO", "KEGG")
PLIST<-list()
for (i in 1:length(DL)) {
  DF<-as.data.frame(DL[[i]])
  DF$Pathway<-factor(DF$Term, levels=unique(DF$Term))
  DF$Module<-factor(DF$Module, levels=unique(DF$Module))
  #Calculating ratio of genes participating (DF$DE) in the complete list of genes (DF$N) in every pathway
  DF$GeneRatio<-DF$DE/DF$N
  
  PLIST[[i]]<-ggplot(DF, aes(x= Module, y= Pathway, fill=GeneRatio))+
    geom_tile()+
    #setting white background
    theme_classic()+
    #Choosing colour and title
    scale_fill_gradientn(colours = paleta)+
    ggtitle(paste(base[i]," pathways"))+
    #x thick marks orientation and vertical dashlines 
    theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.8)),
          panel.border = element_rect(fill=NA),
          axis.title = element_blank(),
          axis.text.y = element_text(size = rel(0.8)),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank())
}
#Lado a lado
cowplot::plot_grid(PLIST[[1]], PLIST[[2]],ncol = 2)
# ============================================================================
# Section 5: DESICION
# ============================================================================
#REMOVE (IF NECCESSARY):
#IGLV3 (Plasmatic B cell)
#BIRC3 (Mature DC)
#CD1C (cDC2)
#Labeling cells with kmeans, dividing in 4 groups (Low, Intermediate 1, Intermediate 2, High)
cell_contamination<-c("IGLV3.1.score","BIRC3.score","CD1C.score")
remove_cells<-list()
for (i in 1:length(cell_contamination)){
  #for each signature, separate cells in 4 groups with kmeans
  cont_sub<-Merg.anchors.scores@meta.data[,cell_contamination[i]]
  if(cell_contamination[i]=="CD1C.score"){
    sep<-kmeans(cont_sub, 5)
  }else{
    sep<-kmeans(cont_sub, 4)
  }
  #find cluster with the higher value
  cl<-names(sep$centers[which.max(sep$centers), ])
  remove_cells[[i]]<-rownames(Merg.anchors.scores@meta.data[sep$cluster==cl,])
}
#adding doublets previously selecte
remove_cells[[4]]<-Doublets
#Adding names
names(remove_cells)<-c("Doublets","Mature DC","cDC2","Doublets")
#Removing cells
#to_remove<-unique(as.character(unlist(remove_cells)))
DF.remove<-plyr::ldply(remove_cells, cbind)

#Renaming doublets
BRCA_Refined <- SetIdent(Annotated, cells = DF.remove[,2], value = DF.remove[,1])
BRCA_Refined$celltypes<-Idents(BRCA_Refined)
#Saving
#saveRDS(BRCA_Refined, file = "D:/R_ScriptsPaper/Objects/WuBRCA_Refined.rds")
# ============================================================================
# Section 6: SECOND TAM RECLUSTERING
# ============================================================================
#List of Seurat Objects with no normalization, clustering or annotation 
my_seurat<-readRDS(file = "D:/Munster/Muenster Workstation/huertaR/Objects/Seurat_V4/Wu_Seurat_V4.rds")
my_seurat@meta.data[,c("celltype_subset","celltype_minor","celltype_major")]<-NULL

#Select only TAMs from BRCA_Refined
BRCA_Refined<-readRDS(file = "D:/R_ScriptsPaper/Objects/WuBRCA_Refined.rds")
TAMs.ID<-WhichCells(BRCA_Refined, idents = c("TAMs","Cycling TAMs"))
#Saved idents
TAMidents<-Idents(BRCA_Refined[,TAMs.ID])

S.Obj<-my_seurat[,TAMs.ID]
#Some samples have less than 200 cells. We have to group them in pairs
samp.size<-sort(table(S.Obj$orig.ident),decreasing=T)
samp.name<-names(samp.size)[samp.size<=233]
#Half the samples
first_half<-samp.name[1:(length(samp.name)/2)]
second_half<-rev(samp.name)[1:(length(samp.name)/2)]
#New name in pairs
fir_sec<-paste(first_half,second_half,sep = "_")
#New sample grouping
S.Obj$TAMsample<-as.character(S.Obj$orig.ident)
for (i in 1:length(first_half)) {
  S.Obj$TAMsample<-ifelse(as.character(S.Obj$TAMsample)%in%c(first_half[i],second_half[i]),
                          fir_sec[i],
                          as.character(S.Obj$TAMsample))
}
#Batch correction
#1-Split by sample
Merg_list<-SplitObject(S.Obj, split.by = "TAMsample")
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

#Adding saved idents
Merg.anchors <- SetIdent(Merg.anchors, cells = names(TAMidents), 
                         value = as.character(TAMidents))
Merg.anchors$celltypes<-Idents(Merg.anchors)
#Saving
#saveRDS(Merg.anchors, file = "D:/R_ScriptsPaper/Def_Objects/WuTAMs.rds")