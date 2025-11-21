# =============================================================================
# Name: Manual Annotation Refinement
# Author: elopez
# Date: Jan-16-2024
# Description: SingleR allows annotation cell by cell. Performed in X script
# Combined with manual anotation
# Designed for Seurat version 4.3.0. Changes are expected for most recent version

# TODO: Use previously curated singleR references (singleR_ReferenceMatrices)
# for cleaner annotation instead of default references
# =============================================================================
library(Seurat)
library(ggplot2)
#Load clustered Seurat object
Annotated<-readRDS(file = "D:/R_Scripts/EssentialObjects/Annotated.rds")
Idents(Annotated)<-Annotated$seurat_clusters
#Adding sufix
Annotated$celltypes<-paste("SnR.",Annotated$celltypes, sep = "")
Annotated$Annotation<-Annotated$celltypes

#Color palette for FeaturePlots
reds<-RColorBrewer::brewer.pal(9,"YlOrRd")
# ============================================================================
# Section 1: Automatic annotation with singleR
# ============================================================================
#Unify counts, normalized data and scale data into a single RNA assay
Annotated[["RNA"]] <- JoinLayers(Annotated[["RNA"]])

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
saveRDS(pred.hesc, file = "/home/huerta/huertaR/Objects/singleR_results.rds")

#Converting to a df
Clust.annot<-as.data.frame(cbind(rownames(pred.hesc),pred.hesc$labels))
#Assigning cell type lables to cells
new.cluster.ids <- Clust.annot$V2
names(new.cluster.ids) <- Clust.annot$V1
#Save data in metadata
Annotated$Annotation<-new.cluster.ids
#
# ============================================================================
# Section 2: MANUAL REFINEMENT OF ANNOTATION: DEGs
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

# ============================================================================
# Section 2.1: REFINEMENT AND CONFIRMATION OF LYMPHOID CELLS
# ============================================================================
#Markers for lymphoid populations
LyMark<-c("CD3E","CD3D", "TRAC",#T CELLS,
          "CD8A","CD8B","GZMA", #CD8
          "IL7R","CD4","CCR7","LEF1", #CD4
          "CTLA4","FOXP3", #T regs
          "NKG7","KLRD1", #NK
          "MS4A1","CD19", #B CELLS
          "IGHG4","IGHG1","CD38","SSR4","JCHAIN","IGHA1" #PLASMATIC B CELL
)
#Potential Lymphoid clusters
Ly.Cl<-unique(as.character(top15$cluster[top15$gene%in%LyMark]))
#Visual confirmation of lymphoid clusters
pA<-DimPlot(Annotated, group.by = "celltypes")
pB<-DimPlot(Annotated, label = T) + NoLegend()
pA+pB
#Expression by cluster
VlnPlot(Annotated, features = LyMark, group.by = "seurat_clusters", 
        assay = "RNA",pt.size = 0)
print(Ly.Cl)

#Lymphoid cells have the following characteristics
Lymph<-subset(Annotated,
              #First lymphoid cells
              subset=Annotation%in%c("SnR.T_cells","SnR.B_cell","SnR.NK_cell")&
                #Second select clusters expressing lymphoid markers
                seurat_clusters%in%c(Ly.Cl,0,1,3,12))

#Visualize expression
p1<-DimPlot(Lymph, group.by = "celltypes")
p2<-DimPlot(Lymph, label = T)
p1+p2
#Corroborate selected genes
FeaturePlot(Lymph, features = "SSR4",slot = "data",cols = reds)

# CONCLUSION ACCORDING TO MARKER EXPRESSION PATTERN:
# NK are correctly annotated
# No plasma B cells subset
# T cell must be refined into plasma B cell and CD4/CD8/T regs

#Defining CD8
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c("8","1")&
                               Annotated$Annotation=="SnR.T_cells",
                             "T CD8",
                             as.character(Annotated$Annotation))
#Defining CD4
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c("0","3","12")&
                               Annotated$Annotation=="SnR.T_cells",
                             "T CD4",
                             as.character(Annotated$Annotation))
#Defining Tregs
Annotated$Annotation<-ifelse(Annotated$seurat_clusters=="11"&
                               Annotated$Annotation=="SnR.T_cells",
                             "T regs",
                             as.character(Annotated$Annotation))

#Defining NK
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(7,10)&
                               Annotated$Annotation=="SnR.NK_cell",
                             "NK cell",
                             as.character(Annotated$Annotation))

#Defining B cells
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(6)&
                               Annotated$Annotation=="SnR.B_cell",
                             "B cell",
                             as.character(Annotated$Annotation))

# ============================================================================
# Section 2.2: REFINEMENT AND CONFIRMATION OF MYELOID CELLS
# ============================================================================
#Markers for myeloid populations
MMark<-c("LYZ", #Myeloid
         "CD68","CD163", "CD14","FCGR3A", #Macrophages y monocytes (TAMs)
         "CPVL","CLEC9A","LGALS2", #cDC1
         "CD1C","FCER1A","SIRPA", #cDC2
         "CCR7","LAMP3","CCL17",#Mature DC
         "CXCR3","LILRA4", #pDC
         "FCGR3B","CXCR2", #Neutrophils
         "TPSAB1","CPA3" # Mast cells
)

#Potential myeloid clusters
Myel.Cl<-unique(as.character(top15$cluster[top15$gene%in%MMark]))
pA+pB
#Visual confirmation
VlnPlot(Annotated, features = MMark, group.by = "seurat_clusters", 
        assay = "RNA",pt.size = 0)
print(Myel.Cl)

#Myeloid cells have the following characteristics
Myeloid<-subset(Annotated,
                #First myeloid cells
                subset=Annotation%in%c("SnR.Macrophage","SnR.Monocyte","SnR.DC",
                                       "SnR.Neutrophils","SnR.Mast_cell")&
                  #Second select seurat clusters of myeloid identity
                  seurat_clusters%in%Myel.Cl[c(-2)])


#Corroborate selected genes
FeaturePlot(Myeloid, features = c("CPVL","CD1C","FCER1A","LAMP3"),slot = "data",cols = reds)

#CONCLUSION ACCORDING TO MARKER EXPRESSION PATTERN
# Merge Macrophages and Monocytes as Mac/Mo
# cDC1 are indistinguishable from Macs
# Only cDC2 and Mature LAMP3 markers are detectable. Comprise in "DC populations" group

#CHANGE CELL TYPES
#Macrophages, monocytes are TAMs except those expressing DC markers
Annotated$Annotation<-ifelse(Annotated$Annotation%in%c("SnR.Macrophage","SnR.Monocyte",
                                                       "SnR.DC"),
                             #&Annotated$seurat_clusters%in%c(5),
                             "TAMs",
                             as.character(Annotated$Annotation))

#Discriminating between blood monocytes and Tissue resident macrophages
Annotated$Annotation<-ifelse(Annotated$Annotation=="TAMs"&grepl("_B",Annotated$Samp_Cond),
       "Blood monocytes",
       ifelse(Annotated$Annotation=="TAMs"&grepl("_N",Annotated$Samp_Cond),
              "Tissue macrophages",
              as.character(Annotated$Annotation)))

#Defining pDC
Annotated$Annotation<-ifelse(Annotated$seurat_clusters==15|
                               Annotated$Annotation=="SnR.pDC",
                             "pDC",
                             as.character(Annotated$Annotation))

#Defining Neutrophils
Annotated$Annotation<-ifelse(Annotated$seurat_clusters==13&
                               Annotated$Annotation=="SnR.Neutrophils",
                             "Neutrophils",
                             as.character(Annotated$Annotation))

#Defining Mast cells
Annotated$Annotation<-ifelse(Annotated$seurat_clusters==14&
                               Annotated$Annotation=="SnR.Mast_cell",
                             "Mast cells",
                             as.character(Annotated$Annotation))

#Defining DC populations
DCs<-WhichCells(Myeloid, expression = CD1C >= 1|FCER1A >= 1|LAMP3>=1)
Annotated$Annotation<-ifelse(colnames(Annotated)%in%DCs&
                               Annotated$seurat_clusters%in%c("5"),
                             "DC populations",
                             as.character(Annotated$Annotation))
# ============================================================================
# Section 2.3: REFINEMENT OF FIBROBLAST,EPITHELIAL,ENDOTHELIAL CELLS ANNOTATION
# ============================================================================
FE_Mark<-c("COL1A1","DCN",#Fibroblast
           "RGS5","ACTA2", #Pericytes
           "EPCAM","KRT18", #Epithelial_cells
           "PECAM1","VWF") #Endothelial_cells
End_Epi_Fib<-subset(Annotated,subset=Annotation%in%c("SnR.Endothelial_cells","SnR.Epithelial_cells",
                                                     "SnR.Fibroblasts","SnR.Smooth_muscle_cells",
                                                     "SnR.Tissue_stem_cells","SnR.Chondrocytes"))

#Visualize expression
pC<-DimPlot(End_Epi_Fib, group.by = "celltypes")+theme(legend.position = "bottom")
pA<-DimPlot(End_Epi_Fib, group.by = "seurat_clusters",label = T)+NoLegend()
pB<-FeaturePlot(End_Epi_Fib, features = FE_Mark,slot = "data",cols = reds)
(pC/pA)|pB

VlnPlot(End_Epi_Fib, features = FE_Mark, group.by = "seurat_clusters", 
        assay = "RNA",pt.size = 0)

#CONCLUSION ACCORDING TO MARKER EXPRESSION PATTERN
# Endothelial, epithelial are correctly annotated.
# All cells from cluster 9 are epithelial
# No pericytes
# Smooth_muscle_cells, Tissue_stem_cells and Chondrocytes must be grouped into fibroblast

#Defining Fibroblasts
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(2,4)&
                               Annotated$Annotation%in%c("SnR.Fibroblasts","SnR.Smooth_muscle_cells",
                                                     "SnR.Tissue_stem_cells","SnR.Chondrocytes"),
                             "Fibroblasts",
                             as.character(Annotated$Annotation))
#Defining Epithelial cells
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(9)&
                               Annotated$Annotation%in%c("SnR.Epithelial_cells"),
                             "Epithelial cells",
                             as.character(Annotated$Annotation))
#Defining Endothelial cells
Annotated$Annotation<-ifelse(Annotated$seurat_clusters%in%c(2,4,5)&
                               Annotated$Annotation%in%c("SnR.Endothelial_cells"),
                             "Endothelial cells",
                             as.character(Annotated$Annotation))

# ============================================================================
# Section 2.4: IDENTIFYING CYCLING CELLS (IF ANY)
# ============================================================================
pCy1<-FeaturePlot(Annotated, features = c("MKI67","CDK1"),cols = reds)
pCy2<-DimPlot(Annotated, group.by = "seurat_clusters", label = T)+NoLegend()
pCy2|pCy1
#CONCLUSION ACCORDING TO MARKER EXPRESSION PATTERN
#No cycling clusters
# ============================================================================
# Section 2.5: RE-ASSIGNING CELLS LEFT
# ============================================================================
#Correcting cells with the annotation from singleR
Rare<-unique(Annotated$Annotation[grepl("SnR.",Annotated$Annotation)])
Rare_IDs<-colnames(Annotated)[Annotated$Annotation%in%Rare]
#Visualize distribution of rare cell types
DimPlot(Annotated, cells.highlight = Rare_IDs)

#------Re-annotate rare cell types with their respective seurat cluster---------
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

#Save object with refined identities
#saveRDS(Annotated, file = "D:/R_ScriptsPaper/Objects/AzziAnnotated.rds")
#
# ============================================================================
# Section 3: Visualization
# ============================================================================
library(ggplot2)
library(cowplot)
#Plot batch-corrected samples
p.type<-DimPlot(Annotated, reduction = "umap",group.by = "celltypes")+
  labs(title = "SingleR")+
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
Disc.Pal<-DiscretePalette(length(unique(Annotated$Annotation)),shuffle = T)

p.annot<-DimPlot(Annotated, reduction = "umap",label = F,
                 group.by = "Annotation",cols = Disc.Pal)+
  labs(title = "Annotation")+
  theme(axis.title = element_text(size = rel(0.6)),
        axis.text = element_text(size = rel(0.6)))

#All 3 plots
p.type+p.cl+p.annot

#Final annotation (DotPlot)
AllMarkers<-c("CD3D","CD3E",
              "CD8A",
              "CD4","IL7R",
              "FOXP3","CTLA4",
              "NKG7","KLRD1",
              "MS4A1","CD19",
              "CD68","CD14","FCGR3A",
              "CLEC9A","CD1C","FCER1A","LAMP3", #DC populations
              "FCGR3B","CXCR2", #Neutrophils
              "TPSAB1","CPA3", # Mast cells
              "CXCR3","LILRA4",
              "COL1A1","DCN", #Fibroblast
              "VWF","PECAM1", #Endothelial cells
              "EPCAM","KRT19" #Epithelial_cells
              )

#Order on Y axis
type.order<-c("T CD8","T CD4","T regs","NK cell","B cell",
          "TAMs","Tissue macrophages","Blood monocytes","DC populations",
          "Neutrophils","Mast cells","pDC",
          "Fibroblasts",
          "Endothelial cells","Epithelial cells")

DotPlot(Annotated, features = AllMarkers, assay = "RNA", cols = "Spectral")+
  scale_y_discrete(limits=rev(type.order))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = rel(0.8)),
        axis.text.y = element_text(size = rel(0.8)),
        axis.title = element_blank())
