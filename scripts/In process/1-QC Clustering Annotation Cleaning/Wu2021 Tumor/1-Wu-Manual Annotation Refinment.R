# =============================================================================
# Name: Annotation
# Author: huerta
# Date: Jan-16-2024
# Description: Corroboration of annotated cell types by Wu, 2021 with singleR.
# SingleR allows annotation cell by cell

# TODO: Use previously curated singleR references (singleR_ReferenceMatrices)
# for cleaner annotation instead of default references
# =============================================================================
library(Seurat)
library(ggplot2)
#Load clustered Seurat object
Annotated<-readRDS(file = "D:/Munster/Muenster Workstation/huertaR/Objects/Seurat_V4/SnRAnnotWu_V4.rds")

#Assign as identities
Idents(Annotated)<-Annotated$seurat_clusters

#Adding sufix
Annotated$celltypes<-paste("SnR.",Annotated$celltypes, sep = "")
Annotated$Annotation<-Annotated$celltypes

#Color palette for FeaturePlots
reds<-RColorBrewer::brewer.pal(9,"YlOrRd")

#All markers from all populations that will be used
AllMarkers<-c("CD3E","CD3D", "TRAC",#T CELLS,
              "CD8A","CD8B","GZMA", #CD8
              "IL7R","CD4","CCR7", #CD4
              "CTLA4","FOXP3", #T regs
              "NKG7","KLRD1", #NK
              "MS4A1","CD19", #B CELLS
              "IGHG4","IGHG1","CD38","SSR4","IGHA1", #PLASMATIC B CELL
              "LYZ", #Myeloid
              "CD68","CD163", "CD14","FCGR3A", #Macrophages y monocytes (TAMs)
              "CPVL","CLEC9A","XCR1", #cDC1
              "CD1C","FCER1A","CLEC10A", #cDC2
              "CCR7","LAMP3","BIRC3",#Mature DC
              "CXCR3","LILRA4","JCHAIN", #pDC
              "FCGR3B","CXCR2", #Neutrophils
              "TPSAB1","CPA3", # Mast cells
              "COL1A1","DCN",#Fibroblast
              "RGS5","ACTA2", #Pericytes
              "EPCAM","KRT18", #Epithelial_cells
              "PECAM1","VWF", #Endothelilal cells
              "TUBB1","PPBP" #Platelets
)
# ============================================================================
# Section 1: Automatic annotation with singleR
# ============================================================================
Annotated<-readRDS(file = "D:/Munster/Muenster Workstation/huertaR/Objects/Seurat_V4/IntegratedWu_V4.rds")
Annotated@meta.data[,c("celltype_subset","celltype_minor","celltype_major")]<-NULL

#Unify counts, normalized data and scale data into a single RNA assay
#Annotated[["RNA"]] <- JoinLayers(Annotated[["RNA"]])

#Extract normalize expression matrix
hESCs<-GetAssayData(Annotated, assay="RNA", layer='data')
#Removing structural genes
bad_patterns <- 'RP[0-9]|^MT|[0-9][0-9][0-9][0-9]|^RPL|^RPS'
a <- rownames(hESCs)[!grepl(bad_patterns, rownames(hESCs))]
hESCs<-hESCs[a,]

#Automatic annotation per cell 
library(SingleR)
#Load reference transcriptome data from cell lines
References<-readRDS(file = "/home/huerta/huertaR/previous work/Objects/singleR_ReferenceMatrices.rds")

pred.hesc <- SingleR(test = hESCs, 
                     #Reference transcriptomes
                     ref = References,
                     #cell by cell, not per cluster
                     method = "single",
                     assay.type.test=1,
                     labels = References$label.main)
#Summary of annotated cells
df.labels<-as.data.frame(table(pred.hesc$pruned.labels))
library(dplyr)
#If frequency is less than 100. label as unknown
unknown<-as.character(df.labels[df.labels$Freq<=100,1])
new.labels<-ifelse(pred.hesc$pruned.labels%in%unknown,"unknown",pred.hesc$pruned.labels)
pred.hesc$labels<-new.labels

#Visualize preliminary results of SingleR
plotScoreHeatmap(pred.hesc)

#Save singleR results
saveRDS(pred.hesc, file = "D:/Munster/Muenster Workstation/huertaR/Objects/Seurat_V4/singleR_results.rds")

#Converting to a df
Clust.annot<-as.data.frame(cbind(rownames(pred.hesc),pred.hesc$labels))
#Assigning cell type lables to cells
new.cluster.ids <- Clust.annot$V2
names(new.cluster.ids) <- Clust.annot$V1
#Save data in metadata
Annotated$celltypes<-new.cluster.ids
saveRDS(Annotated,file = "D:/Munster/Muenster Workstation/huertaR/Objects/Seurat_V4/SnRAnnotWu_V4.rds")
#
# ============================================================================
# Section 2: MANUAL REFINEMENT OF ANNOTATION AND MEDIANS: DEGs
# ============================================================================
#Calculate marker genes of each cluster
Cl.Mk<-FindAllMarkers(Annotated,
                      #Calculate on normalize data
                      assay = "RNA", slot = "data", 
                      #Minimal fold change and expression %
                      logfc.threshold = 1, min.pct = 0.3, 
                      only.pos=T)

#Top 15 genes per cluster
library(dplyr)
Cl.Mk %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 2) %>%
  slice_head(n = 15) %>%
  ungroup() -> top15

#CALCULATING MEDIANS FOR EACH CLUSTER
S.Obj<-SplitObject(Annotated, split.by = "seurat_clusters")
med.list<-list()
for (i in 1:length(S.Obj)) {
  med.matrix<-GetAssayData(S.Obj[[i]],assay = "RNA",slot = "data")
  #Genes present in the matrix
  SubMark<-AllMarkers[AllMarkers%in%rownames(med.matrix)]
  #Calculate medians
  #med.list[[i]] <- apply(med.matrix[SubMark,], 1, median)
  #Calculate upper quantile
  med.list[[i]] <- apply(med.matrix[SubMark, ], 1, function(x) quantile(x, probs = 0.75))
  print(paste("Succesfull loop",i, sep = " "))
}
names(med.list)<-names(S.Obj)
#Storign results in df
med.df<-plyr::ldply(med.list, rbind)
# ============================================================================
# Section 2.1: REFINEMENT AND CONFIRMATION OF LYMPHOID CELLS
# ============================================================================
#Markers for lymphoid populations
LyMark<-c("CD3E","CD3D", "TRAC",#T CELLS,
          "CD8A","CD8B","GZMA", #CD8
          "IL7R","CD4","CCR7", #CD4
          "CTLA4","FOXP3", #T regs
          "NKG7","KLRD1", #NK
          "MS4A1","CD19", #B CELLS
          "IGHG4","IGHG1","CD38","SSR4","IGHA1" #PLASMATIC B CELL
)
#Potential Lymphoid clusters
Ly.Cl<-unique(as.character(top15$cluster[top15$gene%in%LyMark]))
#Visual confirmation of lymphoid clusters
pA<-DimPlot(Annotated, group.by = "celltypes")
pB<-DimPlot(Annotated, label = T) + NoLegend()
pA+pB
#Expression by cluster
#Ordering clusters
Ly_test<-LyMark[LyMark%in%rownames(med.matrix)]
lista_violines<-list()
for (i in 1:length(Ly_test)) {
  
  VlnOrder<-med.df[,1][order(med.df[,Ly_test[i]], decreasing = TRUE)]
  
  lista_violines[[i]]<-VlnPlot(Annotated, features = Ly_test[i], pt.size=0)+
    #Orden eje x
    scale_x_discrete(limits=VlnOrder)+
    labs(title = Ly_test[i])+
    #El texto orizontal que etiqueta cada curva
    theme(
      axis.text.y = element_text(size = rel(0.5)),
      axis.text.x = element_text(size = rel(0.5)),
      #El texto horizontal
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(size = rel(1), family = "serif")
    )+
    NoLegend()
  print(paste("Succesfull loop",i, sep = " "))
}
#Arreglo violines
gridExtra::grid.arrange(grobs=lista_violines, ncol=5)
print(Ly.Cl)

#Lymphoid cells have the following characteristics
Lymph<-subset(Annotated,
              #First lymphoid cells
              subset=Annotation%in%c("SnR.T_cells","SnR.B_cell","SnR.NK_cell")&
                #Second select clusters expressing lymphoid markers
                seurat_clusters%in%c(Ly.Cl,18))

#Visualize expression
p1<-DimPlot(Lymph, group.by = "celltypes")
p2<-DimPlot(Lymph, label = T)
p1+p2
#Corroborate selected genes
FeaturePlot(Lymph, features = LyMark[7:length(LyMark)],slot = "data",cols = reds)

# CONCLUSION ACCORDING TO MARKER EXPRESSION PATTERN:
# B cells are correctly annotated
# No plasma B cells subset
# T cell must be refined into NK CD4/CD8/T regs by marker expression. 
# All lymphoid cells are in cluster 1

#Defining CD8
CD8<-WhichCells(Lymph, expression = CD8A >= 1|CD8B >= 1)
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(4,25)&
                               Annotated$Annotation=="SnR.T_cells"&
                             colnames(Annotated)%in%CD8,
                             "T CD8",
                             as.character(Annotated$Annotation))
#Defining CD4
CD4<-WhichCells(Lymph, expression = IL7R >= 1)
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(4,5,12)&
                               Annotated$Annotation=="SnR.T_cells"&
                             colnames(Annotated)%in%CD4,
                             "T CD4",
                             as.character(Annotated$Annotation))
#Defining Tregs
tregs<-WhichCells(Lymph, expression = CTLA4 >= 1|FOXP3 >= 1)
Annotated$Annotation<-ifelse(Annotated$seurat_clusters==12&
                               Annotated$Annotation=="SnR.T_cells"&
                               colnames(Annotated)%in%tregs,
                             "T regs",
                             as.character(Annotated$Annotation))

#Defining NK
nk<-WhichCells(Lymph, expression = NKG7 >= 1)
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(4,25)&
                               Annotated$Annotation=="SnR.NK_cell"&
                             colnames(Annotated)%in%nk,
                             "NK cell",
                             as.character(Annotated$Annotation))

#Defining B cells
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(23)&
                               Annotated$Annotation=="SnR.B_cell",
                             "B cell",
                             as.character(Annotated$Annotation))
#Defining plasmatic B cells
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(20),
                             "Plasmatic B cell",
                             as.character(Annotated$Annotation))
# ============================================================================
# Section 2.2: REFINEMENT AND CONFIRMATION OF MYELOID CELLS
# ============================================================================
#Markers for myeloid populations
MMark<-c("LYZ", #Myeloid
         "CD68","CD163", "CD14","FCGR3A", #Macrophages y monocytes (TAMs)
         "CPVL","CLEC9A","XCR1", #cDC1
         "CD1C","FCER1A","CLEC10A", #cDC2
         "CCR7","LAMP3","BIRC3",#Mature DC
         "CXCR3","LILRA4","JCHAIN", #pDC
         "FCGR3B","CXCR2", #Neutrophils
         "TPSAB1","CPA3" # Mast cells
)

#Potential myeloid clusters
Myel.Cl<-unique(as.character(top15$cluster[top15$gene%in%MMark]))
pA<-DimPlot(Annotated, group.by = "Annotation")
pA+pB

#Expression by cluster
#Ordering clusters
MM_test<-MMark[MMark%in%rownames(med.matrix)]
lista_violines<-list()
for (i in 1:length(MM_test)) {
  
  VlnOrder<-med.df[,1][order(med.df[,MM_test[i]], decreasing = TRUE)]
  
  lista_violines[[i]]<-VlnPlot(Annotated, features = MM_test[i], pt.size=0)+
    #Orden eje x
    scale_x_discrete(limits=VlnOrder)+
    labs(title = MM_test[i])+
    #El texto orizontal que etiqueta cada curva
    theme(
      axis.text.y = element_text(size = rel(0.5)),
      axis.text.x = element_text(size = rel(0.5)),
      #El texto horizontal
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(size = rel(1), family = "serif")
    )+
    NoLegend()
  print(paste("Succesfull loop",i, sep = " "))
}
#Arreglo violines
gridExtra::grid.arrange(grobs=lista_violines, ncol=5)
print(Myel.Cl)

#Myeloid cells have the following characteristics
Myeloid<-subset(Annotated,
                #First myeloid cells
                subset=Annotation%in%c("SnR.Macrophage","SnR.Monocyte","SnR.DC",
                                       "SnR.Neutrophils")&
                  #Second select seurat clusters of myeloid identity
                  seurat_clusters%in%c(Myel.Cl,19,24))


#Corroborate selected genes
FeaturePlot(Myeloid, features = c("CD68","CD163", "CD14","FCGR3A", #Macrophages y monocytes (TAMs)
                                  "CPVL","CLEC9A","XCR1", #cDC1
                                  "CD1C","FCER1A","CLEC10A", #cDC2
                                  "CCR7","LAMP3","BIRC3"#Mature DC
                                  ),slot = "data",cols = reds)

# CONCLUSION ACCORDING TO MARKER EXPRESSION PATTERN

#CHANGE CELL TYPES
#Defining TAMs
Annotated$Annotation<-ifelse(Annotated$Annotation%in%c("SnR.Macrophage","SnR.Monocyte",
                                                       "SnR.Neutrophils","SnR.DC")&
                               Annotated$seurat_clusters%in%c(7,10,15,19,24),
                             "TAMs",
                             as.character(Annotated$Annotation))
#Defining pDC
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(26),
                             "pDC",
                             as.character(Annotated$Annotation))

#Defining DC populations (Remaining cells) 
cDC1<-WhichCells(Myeloid, expression = CLEC9A >= 1)
Annotated$Annotation<-ifelse(colnames(Annotated)%in%cDC1&
                               colnames(Annotated)%in%colnames(Myeloid),
                             "cDC1",
                             as.character(Annotated$Annotation))
#cDC2
cDC2<-WhichCells(Myeloid, expression = CD1C>= 2|FCER1A>= 2|CLEC10A>= 2.5)
Annotated$Annotation<-ifelse(colnames(Annotated)%in%cDC2&
                               colnames(Annotated)%in%colnames(Myeloid),
                             "cDC2",
                             as.character(Annotated$Annotation))
#Mature DC
MDC<-WhichCells(Myeloid, expression = CCR7>= 2.8|BIRC3>= 2.8)
Annotated$Annotation<-ifelse(colnames(Annotated)%in%MDC&
                               colnames(Annotated)%in%colnames(Myeloid),
                             "Mature DC",
                             as.character(Annotated$Annotation))
# ============================================================================
# Section 2.3: REFINEMENT OF FIBROBLAST,EPITHELIAL,ENDOTHELIAL CELLS ANNOTATION
# ============================================================================
FE_Mark<-c("COL1A1","DCN",#Fibroblast
           "RGS5","ACTA2", #Pericytes
           "EPCAM","KRT18", #Epithelial_cells
           "PECAM1","VWF", #Endothelial_cells
           "TUBB1","PPBP") #Platelet
FE.Cl<-unique(as.character(top15$cluster[top15$gene%in%FE_Mark]))
#Ordering clusters
FE_test<-FE_Mark[FE_Mark%in%rownames(med.matrix)]
lista_violines<-list()
for (i in 1:length(FE_test)) {
  
  VlnOrder<-med.df[,1][order(med.df[,FE_test[i]], decreasing = TRUE)]
  
  lista_violines[[i]]<-VlnPlot(Annotated, features = FE_test[i], pt.size=0)+
    #Orden eje x
    scale_x_discrete(limits=VlnOrder)+
    labs(title = FE_test[i])+
    #El texto orizontal que etiqueta cada curva
    theme(
      axis.text.y = element_text(size = rel(0.5)),
      axis.text.x = element_text(size = rel(0.5)),
      #El texto horizontal
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      plot.title = element_text(size = rel(1), family = "serif")
    )+
    NoLegend()
  print(paste("Succesfull loop",i, sep = " "))
}
#Arreglo violines
gridExtra::grid.arrange(grobs=lista_violines, ncol=3)
print(FE.Cl)
pA<-DimPlot(Annotated, group.by = "Annotation")
pA+pB
#
End_Epi_Fib<-subset(Annotated,subset=seurat_clusters%in%c(FE.Cl,0:3,8,18,16,22,28,14))

#Visualize expression
pc<-DimPlot(End_Epi_Fib, group.by = "celltypes")+theme(legend.position = "bottom")
pa<-DimPlot(End_Epi_Fib, group.by = "seurat_clusters",label = T)+NoLegend()
pb<-FeaturePlot(End_Epi_Fib, features = FE_Mark,slot = "data",cols = reds)
(pc/pa)|pb

#Defining Fibroblast
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(9,11,13)&
                               Annotated$Annotation%in%"SnR.Fibroblasts",
                             "Fibroblasts",
                             as.character(Annotated$Annotation))

#Defining Pericytes
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(16,17,21)&
                               Annotated$Annotation%in%"SnR.Fibroblasts",
                             "Pericytes",
                             as.character(Annotated$Annotation))
#Defining Endothelial cells
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(6,14,22,27,28)&
                               Annotated$Annotation%in%"SnR.Endothelial_cells",
                             "Endothelial cells",
                             as.character(Annotated$Annotation))
#Defining Epithelial cells
Annotated$Annotation<-ifelse(Annotated$Annotation%in%"SnR.Epithelial_cells",
                             "Epithelial cells",
                             as.character(Annotated$Annotation))
# ============================================================================
# Section 2.4: IDENTIFYING CYCLING CELLS (IF ANY)
# ============================================================================
pCy1<-FeaturePlot(Annotated, features = c("MKI67","CDK1"),cols = reds)
pCy2<-DimPlot(Annotated, group.by = "seurat_clusters", label = T)+NoLegend()
pCy2|pCy1
#CONCLUSION ACCORDING TO MARKER EXPRESSION PATTERN
#clusters 25 and 18 are cycling
#Adding cycling prefix
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(18,24,25)&
                               Annotated$Annotation%in%c("TAMs","T CD8"),
                             paste("Cycling",as.character(Annotated$Annotation),sep = " "),
                             as.character(Annotated$Annotation))
# ============================================================================
# Section 2.5: RE-ASSIGNING CELLS LEFT
# ============================================================================
#Correcting cells with the annotation from singleR
Rare<-unique(Annotated$Annotation[grepl("SnR.",Annotated$Annotation)])
Rare_IDs<-colnames(Annotated)[Annotated$Annotation%in%Rare]
#Visualize distribution of rare cell types
DimPlot(Annotated, cells.highlight = Rare_IDs)

#Re-annotate rare cell types with their respective seurat cluster
Mixture<-ifelse(
  #Cells identified as rare
  Annotated$Annotation%in%Rare|
    #and cells with no annotation (if any)
    is.na(Annotated$Annotation),
  #re-annotate with numeric seurat cluster (provisionally)
  as.character(Annotated$seurat_clusters),
  #otherwise, do not change identity
  as.character(Annotated$Annotation)
)

#Assign cell type that constitute the majority of the respective seurat cluster
cl_names<-levels(Annotated$seurat_clusters)
max.abund<-vector()
for (j in 1:length(cl_names)) {
  #Calculate abundance of each cell type in every seurat cluster
  abundances<-table(Annotated$Annotation[Annotated$seurat_clusters==cl_names[j]])
  #Cell type with the highest abundance in the cluster
  max.abund[j]<-names(abundances)[abundances==max(abundances)]
}
#Subsitution map
names(max.abund)<-cl_names
#Substitute values
Annotated$Annotation<-dplyr::recode(Mixture, !!!max.abund)
Idents(Annotated)<-Annotated$Annotation
#Removing cols
Annotated$celltypes<-Annotated$Annotation
Annotated$Annotation<-NULL

#Save object with refined identities
#saveRDS(Annotated, file = "D:/R_ScriptsPaper/Objects/WuBRCA_Annotated.rds")
#
# ============================================================================
# Section 3: Visualization
# ============================================================================
library(cowplot)
#Plot batch-corrected samples
p.type<-DimPlot(Annotated, reduction = "umap",group.by = "orig.ident")+
  labs(title = "Merged samples")+
  theme(axis.title = element_text(size = rel(0.6)),
        axis.text = element_text(size = rel(0.6)),
        legend.position = "bottom")

#Plot clustering
p.cl<-DimPlot(Annotated, reduction = "umap",label = T,group.by = "seurat_clusters")+
  labs(title = "Clustering")+
  theme(axis.title = element_text(size = rel(0.6)),
        axis.text = element_text(size = rel(0.6)))+ 
  NoLegend()

# Manually refined annotation
# Use palette of discrette colors
Disc.Pal<-DiscretePalette(length(unique(Annotated$celltypes)),shuffle = T)

p.annot<-DimPlot(Annotated, reduction = "umap",label = F,
                 group.by = "celltypes",cols = Disc.Pal)+
  labs(title = "Annotation")+
  theme(axis.title = element_text(size = rel(0.6)),
        axis.text = element_text(size = rel(0.6)))

#All 3 plots
p.type+p.cl+p.annot

#Final annotation (DotPlot)
AllPairMarkers<-c("CD3D","CD3E",
                  "CD8A", #T CD8
                  "CD4","IL7R", #T CD4
                  "FOXP3","CTLA4", #Tregs
                  "NKG7","KLRD1", # NK cell
                  "MS4A1","CD19", # B cell
                  "IGHG4","IGHG1", # Plasmatic B cell
                  "CD68","CD14","FCGR3A", #TAMs
                  "CLEC9A","XCR1", #cDC1
                  "CD1C","CLEC10A", #cDC2
                  "CCR7","LAMP3",#Mature DC
                  "LILRA4","CXCR3", #pDC
                  "COL1A1","DCN", #Fibroblast
                  "RGS5","ACTA2", #Pericytes
                  "VWF","PECAM1", #Endothelial cells
                  "EPCAM","KRT19", #Epithelial cells
                  "MKI67","CDK1", #Cycling
                  "FCGR3B","CXCR2", #Neutrophils
                  "TPSAB1","CPA3" # Mast cells
                  )

#Order on Y axis
type.order<-c("T CD8","T CD4","T regs","NK cell","B cell","Plasmatic B cell",
              "TAMs","cDC1","cDC2","Mature DC","pDC",
              "Fibroblasts","Pericytes","Endothelial cells", "Epithelial cells",
              "Cycling T CD8", "Cycling TAMs")

DotPlot(Annotated, features = AllPairMarkers, assay = "RNA", cols = "Spectral")+
  scale_y_discrete(limits=rev(type.order))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = rel(0.8)),
        axis.text.y = element_text(size = rel(0.8)),
        axis.title = element_blank())
