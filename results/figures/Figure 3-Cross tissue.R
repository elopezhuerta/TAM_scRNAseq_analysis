# ============================================================================
# Name: Fig3-Cross tissue
# Author: elopez
# Date: Jun-30-2024
# Description: 
# TODO: 
# ============================================================================
library(Seurat)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(patchwork)

# Object with phenotypes assigned in all tissue types
Tissue_Pheno <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Tissue_Pheno.rds")

# Grouping cycling TAMs into TAMs
Tissue_Pheno$celltypes<-replace(as.character(Tissue_Pheno$celltypes), 
                           as.character(Tissue_Pheno$celltypes)=="Cycling TAMs", "TAMs")

# 
colnames(Tissue_Pheno@meta.data)<-stringr::str_replace_all(colnames(Tissue_Pheno@meta.data), "_score"," score")
# Extracting scores
scores <- colnames(Tissue_Pheno@meta.data)[grepl("score",colnames(Tissue_Pheno@meta.data))]
single.modules <-stringr::str_remove_all(scores, ".score")

# Colos of tissue types
pal_tis<-c("#FF0066","#CC0033","#33CCFF","#0033FF","black")
#
# ============================================================================
# Figure 3A: UMAP TissueType
# ============================================================================
# Changing name
Tissue_Pheno$TissueType <- stringr::str_replace_all(Tissue_Pheno$TissueType,"NT-tissue","NTMT")

# Set up the order of tissue types
tissue_types <- c("Blood: Healthy", "Blood: BRCA", "NTMT: Healthy", "NTMT: BRCA", "Tumor")

#Fixing order
Tissue_Pheno$TissueType<-factor(Tissue_Pheno$TissueType,
                           levels = tissue_types)

# Initialize an empty list to store the plots
highlight_list <- list()

# Plots highlighting tissue of origin
for (i in seq_along(tissue_types)) {
  tissue <- tissue_types[i]
  highlight_color <- pal_tis[i]
  
  highlight_list[[i]] <- DimPlot(
    Tissue_Pheno, 
    group.by = "TissueType", 
    cols = highlight_color,
    cells.highlight = colnames(Tissue_Pheno)[Tissue_Pheno$TissueType == tissue],
    sizes.highlight = 0.2,
    pt.size = 0.2
  ) +
    scale_color_manual(
      values = c("lightgrey", highlight_color)
    ) +
    labs(title = tissue)+
    NoLegend()+
    theme(plot.title = element_text(hjust = 0.5, size = 13, face = "plain"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank()
          )
}
# Arrenging plots
UMAP.High <- ggpubr::ggarrange(plotlist = highlight_list, ncol = 1, nrow = 5)
# ============================================================================
# Figure 3B: SCORE PLOTS
# ============================================================================
# Palette color
spect.pal <-rev(RColorBrewer::brewer.pal(11, "Spectral"))

UMAP_Score <- FeaturePlot(Tissue_Pheno, features = scores, max.cutoff = "q99", 
                          cols=spect.pal, pt.size = 0.8, ncol = 7)&
  xlab(NULL)& ylab(NULL)&
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "plain"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

UMAP_Score <- UMAP_Score + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
# ============================================================================
# Figure 3C: VLNPLOTS: ENRICHMENT IN EACH TISSUE
# ============================================================================
# Function to indicate median point in violin
median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}

Vln_scores <- VlnPlot(Tissue_Pheno, features = scores, group.by = "TissueType", 
                      pt.size = 0, ncol = 4, cols = pal_tis)&
  ylab("ES") &
  # More limit to add significant comparisons
  ylim(-1.1,1.3) &
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "plain"),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_blank()
        )&
  # Add point at median
  stat_summary(fun = median.stat, geom='point', size = 1, colour = "white")

# Add a common legend
Vln_scores <- Vln_scores + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12), # Smaller text size for legend
        # Adjust size of the boxes in the label
        legend.key.size = unit(0.4, "cm"))
#
# ============================================================================
# Figure 3D: TRAJECTORY SLINGSHOT
# ============================================================================
# Import seurat with Main Phenotype
MainPhenotype <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/MainPhenotype.rds")
# Extract UMAP dimensions
dimred <- MainPhenotype@reductions$umap@cell.embeddings

#Aesthetic
#Specifying root Monocytes_FCN1
Visualization <- ifelse(MainPhenotype$Tissue=="Blood"&
                          MainPhenotype$MainPheno=="FCN1",
                        "Monocytes FCN1",
                        as.character(MainPhenotype$MainPheno))

# Add labels for visualization
umap.df <-data.frame(dimred, "Phenotype" = Visualization)

#Dicrete palette colors
#custom_palette<-DiscretePalette(length(unique(umap.df$TAM_subset)),shuffle = T)

custom_palette <- c("#FFE100", #APOE
                    "#0075DC", #CXCL9
                    "#2BCE48", #FCN1
                    "#FF0010", #FOLR2
                    "#993F00", #FOSB
                    "#FF5005", #IL1B
                    "#F0A0FF", #ISG15
                    "#33CCFF", #Mono_FCN1
                    "darkgrey" #Unassigned
)

# UMAP plot displaying Main phenotype
pre_traject <- ggplot(umap.df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(fill = Phenotype), shape = 21, size = 1, col="black") +
  # Set the custom combined palette
  scale_fill_manual(values = c(custom_palette)) +
  # Increase the size of the legend circles
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5, 
      override.aes = list(size = 3), 
      #Legend columns
      ncol = 2)) +
  theme_classic() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 10)
  )

# ADDING TRAJECTORY LINES
library(slingshot)
library(ggrepel)
#Import object with Trajectory results
sds <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/TissuesSling.rds")
# Extract Linages as DF
mst <-slingMST(sds, as.df = TRUE)

# Extract the minimum and maximum order for each lineage. 
# In order to know the start and end in each lineage
library(dplyr)
min_max_data <- mst %>%
  group_by(Lineage) %>%
  summarise(min_order = min(Order, na.rm = TRUE), 
            max_order = max(Order, na.rm = TRUE))

# Merge the minimum and maximum order information back to the original data
mst <- merge(mst, min_max_data, by = "Lineage")

# Identify whether a point is the "Start" or "End"
mst$point_type <- ifelse(mst$Order == mst$min_order, 
                         "Start", 
                         ifelse(mst$Order == mst$max_order, "End", "Intermediate"))

# Add red points for "Start" and green points for "End"
Traject_plot <- pre_traject +
  geom_point(data = mst, aes(color = point_type), size = 3) +
  geom_path(data = mst %>% arrange(Order), aes(group = Lineage),
            color = "#FDDDE6",linewidth = 1) +
  scale_color_manual(values = c("Start" = "blue", "End" = "red", "Intermediate" = "#FDDDE6")) +
  # Change title of second legend
  labs(color = "Place in trajectory") +
  guides(
    color = guide_legend(
      title.position = "top", 
      title.hjust = 0.5,
      override.aes = list(size = 3)
    )
  ) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 10)
  )
#
# ============================================================================
# ARRANGE ALL PLOTS
# ============================================================================
DE <-cowplot::plot_grid(Traject_plot, NULL,ncol = 2,
                       rel_widths = c(1,0.7))

BCDE <-cowplot::plot_grid(UMAP_Score,Vln_scores, DE,nrow = 3, 
                         rel_heights = c(0.5,1,0.8)
                         )

cowplot::plot_grid(UMAP.High,BCDE,ncol = 2,
                   #labels = c("a"), label_size = 20,
                   rel_widths = c(0.2,1)
                   )
#