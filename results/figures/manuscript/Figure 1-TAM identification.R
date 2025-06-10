# ============================================================================
# Name: Figure 1
# Author: elopez
# Date: Oct-08-2024
# Description: 
# TODO: 
# ============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)
#Importar objeto BRCA_Refined
Annotated <- readRDS(file = "D:/R_ScriptsPaper/Objects/WuBRCA_Refined.rds")
#Removing doublets for visualization
Annotated <- subset(Annotated, idents = "Doublets", invert=T)

# reds<-RColorBrewer::brewer.pal(9, "YlOrRd")
spectral.pal <- RColorBrewer::brewer.pal(11, "Spectral")
red.pal <- RColorBrewer::brewer.pal(9, "YlOrRd")
#
# ============================================================================
# FIGURE 1B-C: CELL TYPES AND MARKERS
# ============================================================================
# Generate palette of discrette colors
#Disc.Pal<-DiscretePalette(length(unique(Annotated$celltypes)),shuffle = T)
# Custom palette inspired in the previous function
Disc.Pal <- c("#4C005C","#426600","#993F00","#191919","#E0FF66","#2BCE48","#FFE100","#94FFB5",
             "#00998F","#0075DC","#003380","#FFA8BB","#FFCC99","#F0A0FF","#9DCC00",
             "#990000","#FF0010")

Annotation_Plot <- DimPlot(Annotated, reduction = "umap",label = F,
                 group.by = "celltypes", cols = Disc.Pal) +
  labs(title = "Wu dataset") +
  theme(
    plot.title = element_text(size = 12, face = "plain"),
    axis.title = element_text(size = rel(0.5)),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    #Space between figure bow and text
    legend.spacing.x =  unit(0.1, 'mm'),
    # Adjust the size of the text legend
    legend.text = element_text(size = rel(0.6))
    )

#-----------------DOTPLOT: MARKERS
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
#Doublets are not included
type.order<-c("T CD8","T CD4","T regs","NK cell","B cell","Plasmatic B cell",
              "TAMs","cDC1","cDC2","Mature DC","pDC",
              "Fibroblasts","Pericytes","Endothelial cells", "Epithelial cells",
              "Cycling T CD8", "Cycling TAMs")

p.Dot <- DotPlot(Annotated, features = rev(AllPairMarkers), assay = "RNA", cols = "Spectral")+
  scale_y_discrete(limits=type.order)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.6)),
        axis.text.y = element_text(size = rel(0.6)),
        axis.title = element_blank(),
        legend.title = element_text(size = rel(0.6)),  # Adjust legend title size
        # Adjust legend text size
        legend.text = element_text(size = rel(0.5)))+
  coord_flip()
# ============================================================================
# FIGURE 1D: TAMs AND MARKERS
# ============================================================================
WuTAMs <-readRDS(file="D:/R_ScriptsPaper/Def_Objects/WuTAMs.rds")

#-----------PLOT BY CLUSTER
WuTAMs$seurat_clusters <-paste0("C",WuTAMs$seurat_clusters)

plot_TAMs <- DimPlot(WuTAMs, group.by = "seurat_clusters", 
                     label.size = 4, label = TRUE, repel = TRUE) +  # Increased label.size for larger text
  NoLegend() +
  labs(title = "8581 TAMs (Wu dataset)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),
    axis.title = element_text(size = rel(0.5)),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
#
#-----------MONO MAC MARKERS
#Expression levels for macrophages monocytes markers
featsM <-c("CD68", "CD163","CD14","FCGR3A","HLA-DQB1","HLA-DPB1")

MacMoMarkers <- FeaturePlot(WuTAMs, features = featsM,pt.size = 0.8,
                            min.cutoff = "q01",
                            max.cutoff = "q99",
                            cols= red.pal) &
  theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "plain"),
    axis.line = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

MacMoMarkers <- MacMoMarkers + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
#
# ============================================================================
# ARRENGING PLOTS
# ============================================================================
DE <- cowplot::plot_grid(plot_TAMs, MacMoMarkers, nrow = 2,
                        #labels = c("d","e"), label_size = 26,
                        rel_heights = c(0.8, 1))

cowplot::plot_grid(Annotation_Plot,p.Dot,DE, ncol = 3,
                   #labels = c("b","c"), label_size = 26,
                   rel_widths = c(1,0.8,0.9))
