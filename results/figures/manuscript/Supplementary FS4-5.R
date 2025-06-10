# ============================================================================
# Name: Figure1-PhenoPaths
# Author: elopez
# Date: Jun-30-2024
# Description: 
# TODO: 
# ============================================================================
library(Seurat)
library(ggplot2)
# Truncated color palette to make the soluette evident.
red.pal <- RColorBrewer::brewer.pal(9, "YlOrRd")[3:9]
# ============================================================================
# SUPP.FIG. S4B-C: NTMT CELLTYPES AND FEATURE PLOT LIST
# ============================================================================
# Import Non tumoral mammary tissue
HNormal <- readRDS(file = "D:/R_ScriptsPaper/Objects/Twigger_HNormal.rds")

#UMAP
Normal_UMAP <- DimPlot(HNormal, pt.size = 0.1, 
                       label = T, label.size = 2.5, repel = T)+
  ggtitle("Non-tumoral tissue: Healthy")+NoLegend()+
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
        axis.title = element_text(size = 6),
        axis.text = element_blank(),
        axis.ticks = element_blank())

#Expression levels for macrophages monocytes markers
featsM <- c("CD68","CD163","CD14","FCGR3A")

FeatPlot_list <- FeaturePlot(HNormal, features = featsM,
                             cols = red.pal, pt.size = 0.1, ncol = 2) &
  theme(
    plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.text = element_text(size = 8))

library(patchwork)
FeatPlot_list <- FeatPlot_list + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# Arrenging
AB <- cowplot::plot_grid(NULL, Normal_UMAP, ncol = 2)
UMAP_NTMT <- cowplot::plot_grid(AB, FeatPlot_list, nrow = 2)

jpeg("D:/R_ScriptsPaper/Figures/Manuscrito/pre_SuppFig4b-c.png", 
     res = 600,width = 14, height = 17, units = "cm")
plot(UMAP_NTMT)
dev.off()
#
# ============================================================================
# SUPP.FIG. S5A-F: UMAP BLOOD AND NEUTROPHIL MARKERS
# ============================================================================
#HEALTHY BLOOD
HBlood <- readRDS(file = "D:/R_ScriptsPaper/Objects/Wang_HBlood.rds")
#UMAP
Sana_UMAP <- DimPlot(HBlood,pt.size = 0.1, 
                     label = T, label.size = 2.5, repel = T)+
  NoLegend()+
  ggtitle("Healthy blood")+
  #Legend on Top
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
        axis.title = element_text(size = 6),
        axis.text = element_blank(),
        axis.ticks = element_blank())

#Expression levels for macrophages/monocytes markers
featsMono<-c("CD14","FCGR3A","CXCR2","FCGR3B")

FeatsMono_list <- FeaturePlot(HBlood, features = featsMono,
                             cols = red.pal, pt.size = 0.1, ncol = 2) &
  theme(
    plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.text = element_text(size = 8))

library(patchwork)
FeatsMono_list <- FeatsMono_list +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

# Arrenging
SangreSana <- cowplot::plot_grid(Sana_UMAP, FeatsMono_list, ncol = 2)


#BRCA BLOOD
BrCBlood <- readRDS(file = "D:/R_ScriptsPaper/Objects/Brechbuhl_BrCBlood.rds")

BRCA_UMAP <- DimPlot(BrCBlood,pt.size = 0.1, 
                   label = T,label.size = 2.5, repel = T)+
  NoLegend()+
  ggtitle("BRCA blood")+
  #Legend on Top
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
        axis.title = element_text(size = 6),
        axis.text = element_blank(),
        axis.ticks = element_blank())

FeatsMono_list_BRCA <- FeaturePlot(BrCBlood, features = featsMono,
                              cols = red.pal, pt.size = 0.1, ncol = 2) &
  theme(
    plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.text = element_text(size = 6))

FeatsMono_list_BRCA <- FeatsMono_list_BRCA +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")

SangreBRCA <- cowplot::plot_grid(BRCA_UMAP, FeatsMono_list_BRCA, ncol = 2)

jpeg("D:/R_ScriptsPaper/Figures/Manuscrito/pre_SuppFigS5.png", 
     res = 600,width = 18, height = 17, units = "cm")

cowplot::plot_grid(SangreSana,SangreBRCA,nrow = 2, 
                   rel_widths = 1,rel_heights = 1)
dev.off()

#