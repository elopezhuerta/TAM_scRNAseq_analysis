# ============================================================================
# Name: SuppFigS7
# Author: elopez
# Date: Sep-27-2024
# Description: Comparing Tissue resident phentoypes vs monocyte derived
# TODO: 
# ============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
#Color palette
reds <- RColorBrewer::brewer.pal(9, "YlOrRd")
single.modules <- c("FCN1","ISG15","CXCL9","IL1B","FOSB","FOLR2","APOE")

# Importing TAMs from Wu dataset
TAM_Modules <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_Pheno.rds")
ES.M1M2.TAMs <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/M1M2_ES_WuTAMs.rds")
TAM_M1M2_scores <- AddMetaData(TAM_Modules, ES.M1M2.TAMs)

Poles <- subset(TAM_M1M2_scores, idents = single.modules)
# Setting order
Idents(Poles) <- factor(as.character(Idents(Poles)),
                       levels = single.modules)

#Importing signature to test
df_M1M2 <- readxl::read_excel("D:/R_ScriptsPaper/Tablas/M1M2.xlsx", 
                              col_names = T, range = "A1:B58")
colnames(df_M1M2) <- c("Marker","gene")

# As list
M1_M2 <- list(
  M1_Curated = df_M1M2$gene[df_M1M2$Marker=="M1 marker"],
  M2_Curated = df_M1M2$gene[df_M1M2$Marker=="M2 marker"]
)
#
# ============================================================================
# SUPP FIGURE S6a: ES IN EACH SUBSET
# ============================================================================
# Divide scores in subsets to calculate averga ES
DFL <- dplyr::group_split(Poles@meta.data[,c("Modules",colnames(ES.M1M2.TAMs))], Modules)
mean_ES <- list()
mod_ES <- vector()
# Average per signature
for (g in 1:length(DFL)){
  d.f<-as.data.frame(DFL[[g]][-1])
  mean_ES[[g]]<-colSums(d.f)/nrow(d.f)
  #Creando vector con el nombre
  mod_ES[g]<-unique(as.matrix(DFL[[g]])[,1])
}
names(mean_ES) <- mod_ES
#List to DF
proms <- plyr::ldply(mean_ES,rbind)


colnames(proms) <- c("Subset","M1 curated score","M2 curated score")
promsDF <- reshape2::melt(proms)

#HEATMAP
library(RColorBrewer)
library(stringr)
paleta <- rev(brewer.pal(11,"Spectral"))

#CURATED
ES_PER_SUBSET <- ggplot(promsDF, aes(x = variable, y = Subset, fill = value)) +  # Swap x and y
  #Heatmap manual
  geom_tile() +
  # Reset the background
  theme_classic() +
  # Set aspect ratio
  #coord_fixed(ratio = 1/0.7) +
  # Text inside tiles
  geom_text(label = round(promsDF$value, digits = 2), size = 3.5, colour = "black") +
  # Colour palette
  scale_fill_gradientn(colours = paleta, limits = c(-1, 1)) +
  # Labels
  #ylab(NULL) + xlab(NULL) +
  labs(fill = "Avg ES") +
  # Custom order on X (Now this is for the Y axis because we swapped them)
  scale_y_discrete(limits = rev(single.modules),
                   #Custom names (differ from Yorder)
                   labels=str_replace(rev(single.modules),"_","-")) +
  # Custom order on Y (Now this is for the X axis)
  scale_x_discrete(limits = c("M1 curated score","M2 curated score")) +
  # Customize theme
  theme(
    axis.text.x = element_text(colour = c("red", "blue"),
                               angle = 30, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    #plot.title = element_text(hjust = 0.5, color = "blue",
     #                         size = rel(1.1)),
    legend.title = element_text(hjust = 0.5, face = "plain", size = 10),
    panel.border = element_rect(fill = NA),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank()
  )
# ============================================================================
# SUPP FIGURE S6b: TME M1 M2 SCORES
# ============================================================================
# Import DF with ES
ES.M1M2.WuBRCA <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/M1M2_ES_WuBRCA.rds")
M1M2.WuBRCA <- readRDS(file = "D:/R_ScriptsPaper/Objects/WuBRCA_Refined.rds")
#Removing doublets for visualization
M1M2.WuBRCA <- subset(M1M2.WuBRCA, idents = "Doublets", invert=T)
# Adding score metadata
M1M2.WuBRCA <- AddMetaData(M1M2.WuBRCA, ES.M1M2.WuBRCA)

# Replacing TAM identities
M1M2.WuBRCA <- SetIdent(M1M2.WuBRCA, cells = names(Poles$Modules), 
                        value = paste0("TAM_",as.character(Poles$Modules)))

# Removing remaining TAMs
M1M2.WuBRCA <- subset(M1M2.WuBRCA, idents = c("TAMs","Cycling TAMs"), invert=T)
# Identities should be the same as celltypes
M1M2.WuBRCA$celltypes <- Idents(M1M2.WuBRCA)

scores <- c("M1_Curated_score","M2_Curated_score")
titulo <- c("M1 curated signature","M2 curated signature")

two_colors <- c("red","blue")
curvas<-function(ii){
  iPlot<-RidgePlot(M1M2.WuBRCA, features = scores[ii])+
    xlab("Enrichment score")+labs(title = titulo[ii])+
    #El texto orizontal que etiqueta cada curva
    theme(axis.text = element_text(size = 8),
          #El texto en vertical
          axis.title.y = element_blank(),
          axis.title.x = element_text(hjust = 0.5, size = 10),
          #El titulo
          plot.title = element_text(size=10, colour = two_colors[ii],
                                    hjust = 0.5, face = "plain"))+
    #Sin la leyenda
    NoLegend()
}
lista_Ridge <- lapply(seq_along(scores), FUN = "curvas")

TMEM1M2_Ridge <- cowplot::plot_grid(plotlist =lista_Ridge, ncol = 2)
#
# ============================================================================
# SUPP FIGURE S7c: INDIVIDUAL EXPRESSION OF M1/M2 MARKERS
# ============================================================================
# Display all M1 M2 markers except:
No_Display <- c(
  # Already displayed on Main figure
  "IL1B","CXCL10","TNF","CD80","CD86","IRF1",
  "CD163","MARCO","FN1","IL10","TGFB1","CD274","RETNLB","NOS2",
  # Low expression
  "IL12B",
  "CCL17","CCL22","CCL24","CD200R1","EGF","CTSS","CTSC","CTSH")

AllMk <- as.character(unlist(M1_M2))[!as.character(unlist(M1_M2))%in%No_Display]

new_legends <- paste0("TAM_", single.modules)
names(new_legends) <- single.modules

color_vector <- c(
  "#F8766D",  # A reddish color for FCN1
  "#D89000",  # An orange color for ISG15
  "#A3A500",  # A yellow-green color for CXCL9
  "#39B600",  # A green color for FOSB
  #"#00BF7D",  # A teal color for HLA_II
  "#00BFC4",  # A light blue color for HSPs
  "#619CFF",  # A blue color for C1QA
  "#F564E3"   # A pink/purple color for APOE
)

# Display median point 
median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}

violin_all_M1M2 <- list()
for (ii in 1:length(AllMk)) {
  # Define marker color (red for M1_Curated, blue otherwise)
  Farbe <- ifelse(AllMk[ii] %in% M1_M2$M1_Curated, "red", "blue")
  
  violin_all_M1M2[[ii]] <- VlnPlot(Poles, features = AllMk[ii], pt.size=0, 
                             assay = "RNA",slot = "data")+
    ylab(AllMk[ii])+
    scale_x_discrete(limits=single.modules)+
    scale_fill_manual(values = color_vector,labels = new_legends) +
    # Force single row legend
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +  
    theme(
      axis.title.y = element_text(size = 10, angle = 90, 
                                  colour = Farbe, vjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 8),
      # Adjusted right and left margins
      plot.margin = unit(c(-0.2, 0.2, -0.2, 0.2), "cm"),  
      legend.text = element_text(size = 10),
      plot.title = element_blank())+
    # Add point at median
    stat_summary(fun = median.stat, geom='point', size = 1, colour = "black")
}
#First part
library(ggpubr)
remaining_M1M2 <- ggarrange(plotlist = violin_all_M1M2,nrow=5, ncol = 7,
                             common.legend = TRUE, legend = "top")
#
# ============================================================================
# JOINING
# ============================================================================
SA_B_D <- cowplot::plot_grid(ES_PER_SUBSET,TMEM1M2_Ridge, NULL,
                           ncol = 3, rel_widths = c(0.6,0.8,0.9))

cowplot::plot_grid(SA_B_D, remaining_M1M2, nrow = 2, rel_heights = c(1,0.8))
#