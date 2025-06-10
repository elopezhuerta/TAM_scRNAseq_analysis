# ============================================================================
# Name: MonovsTissue
# Author: elopez
# Date: Sep-26-2024
# Description: Comparing Tissue resident phentoypes vs monocyte derived
# TODO: 
# ============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
#Color palette
reds <- RColorBrewer::brewer.pal(9, "YlOrRd")
single.modules <- c("FCN1","ISG15","CXCL9","IL1B","FOSB","FOLR2","APOE")

#Importing signature to test
df_M1M2 <- readxl::read_excel("D:/R_ScriptsPaper/Tablas/M1M2.xlsx", col_names = T, range = "A1:B58")
colnames(df_M1M2) <- c("Marker","gene")

# As list
M1_M2 <- list(
  M1_Curated = df_M1M2$gene[df_M1M2$Marker=="M1 marker"],
  M2_Curated = df_M1M2$gene[df_M1M2$Marker=="M2 marker"]
)

# Importing TAMs from Wu dataset
TAM_Modules <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_Pheno.rds")
ES.M1M2.TAMs <- readRDS(file = "D:/R_ScriptsPaper/Objects/Temporal/ES_M1M2_WuTAMs.rds")
# Adjust colnames
colnames(ES.M1M2.TAMs) <- c("Curated_M1_score", "Curated_M2_score")
TAM_M1M2_scores <- AddMetaData(TAM_Modules, ES.M1M2.TAMs)

Poles <- subset(TAM_M1M2_scores, idents = single.modules)
# Setting order
Idents(Poles) <- factor(as.character(Idents(Poles)), levels = single.modules)
#
# ============================================================================
# FIGURE 4a: MONOCYTES vs TR-Mac
# ============================================================================
# Mono like populations
pre_Monos <- colnames(TAM_Modules)[grepl("FCN1|ISG15|CXCL9",TAM_Modules$Modules)]
# Tissue resident like populations
pre_Res <- colnames(TAM_Modules)[grepl("FOLR2|APOE",TAM_Modules$Modules)]
# Shared states
pre_Shared <- colnames(TAM_Modules)[grepl("FOSB|IL1B",TAM_Modules$Modules)]

# Find elements unique to pre_Monos
Monos <- setdiff(pre_Monos, pre_Res)
# Find elements unique to pre_Res
Residentes <- setdiff(pre_Res, pre_Monos)
# Find element unique to Shared states
Shared_States <- setdiff(pre_Shared,c(Monos,Residentes)) 

# Setting idents for plotting
NewIdents <- c(rep("Monocyte-like (FCN1, ISG15, CXCL9)",length(Monos)),
             rep("Tissue macrophage-like (FOLR2, APOE)",length(Residentes)),
             rep("Shared states (FOSB, IL1B)", length(Shared_States)))

MonocitResident <- SetIdent(TAM_Modules, cells = c(Monos,Residentes, Shared_States), 
                            value = NewIdents)

# Removing other phenotypes for cleaner visualization
MonocitResident<-subset(MonocitResident, 
                        idents = c("Monocyte-like (FCN1, ISG15, CXCL9)",
                                   "Tissue macrophage-like (FOLR2, APOE)",
                                   "Shared states (FOSB, IL1B)"))

Plot_MonRes <- DimPlot(MonocitResident)+
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        axis.title = element_text(size = 9),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  guides(color=guide_legend(ncol=1, byrow=TRUE, 
                            #Circle size
                            override.aes = list(size=2.5)))
#
# ============================================================================
# FIGURE 4b: INFLAMATORY vs REPAIR
# ============================================================================
UNI <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Leading_edges.rds")

# Tissue repair
repair_patterns <- "WOUNDING|ENDOCYTOSIS|DIFFERENTIATION|COAGULATION|LYSOSOME"
Repair <- unique(unlist(UNI[grep(repair_patterns,names(UNI))]))

# Interferon high
#IFN_Res <- unique(unlist(UNI[grep("VIRAL|INTERFERON",names(UNI))]))
IFN_Alfa <- unique(unlist(UNI[grep("VIRAL|ALPHA|TYPE_I",names(UNI))]))
IFN_Gamma <- unique(unlist(UNI[grep("GAMMA|TYPE_II",names(UNI))]))
# Immune activation
Immune <- unique(unlist(UNI[grep("HUMORAL|KAPPAB|INFLAM",names(UNI))]))


Resum <- list(Immune, IFN_Alfa,IFN_Gamma, Repair)
names(Resum) <- c("Innate_immunity","IFN_I","IFN_II","Tissue_repair")

#Enrichment of customized signatures
TAM.Infl_Repar <- AddModuleScore(TAM_Modules, features = Resum, 
                               assay = "RNA", name = names(Resum))
#Removing names from score names
colnames(TAM.Infl_Repar@meta.data)<-ifelse(grepl("Innat|IFN|repair",colnames(TAM.Infl_Repar@meta.data)),
                                           stringr::str_remove_all(colnames(TAM.Infl_Repar@meta.data),"[0-9]"),
                                           colnames(TAM.Infl_Repar@meta.data))
plot.list <- list()
for (ii in seq_along(Resum)) {
  # Modify the title by replacing "_" with a space
  plot_title <- gsub("_", " ", names(Resum)[ii])
  # Create the plot
  plot.list[[ii]] <- FeaturePlot(TAM.Infl_Repar, features = names(Resum)[ii], 
                                 min.cutoff = "q01",
                                 cols = reds) +
    labs(title = plot_title) +
    NoLegend() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 13, face = "plain")
    )
}
plot_Inf_Repair <- ggpubr::ggarrange(plotlist = plot.list, ncol=2, nrow = 2,
                                   common.legend = TRUE, legend = "right")
#
# ============================================================================
# FIGURE 4c: MONO AND TISSUE MARKERS 
# ============================================================================
# RNA must be default assay.AverageExpression is applied to default assay
#DefaultAssay(Poles)<-"RNA"

TissMonoMarkers <- c("S100A8","CCR2","CD14","FCGR3A",
                     "MAF","MERTK","FOLR2","LYVE1")

new_xaxis <- paste0("TAM_", single.modules) 
names(new_xaxis) <- single.modules


# Were to store violin plots
vln_list <- list()
for (ii in seq_along(TissMonoMarkers)) {
  # Set dot size conditionally
  dots <- ifelse(TissMonoMarkers[ii] == "CCR2", 0.1, 0) 
  
  vln_list[[ii]]<-VlnPlot(Poles, features = TissMonoMarkers[ii], pt.size=dots,
                          assay = "RNA",slot = "data")+
    NoLegend() +
    labs(title = TissMonoMarkers[ii])+
    ylab("Norm. expr.")+
    # Change x-axis labels
    scale_x_discrete(labels = new_xaxis) +
    theme(axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8),
          axis.title.y = element_text(size = 10),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 13, face = "plain", hjust = 0.5),
          # Reduce space between plots
          plot.margin=unit(c(-0.5,-0.5,-0.5,0.5), "cm"))
  
  #Loop tracking
  print(paste("Loop",ii, sep = " "))
}
# Plot array
TissMono_Plot <- ggpubr::ggarrange(plotlist = vln_list, ncol=4, nrow = 2)
#
# ============================================================================
# FIGURE 4d: M1/M2 SCORES
# ============================================================================
M1M2_scores <- colnames(ES.M1M2.TAMs)

#-------------------SCORES ON UMAP
# Plots for M1 M2 ES
plot.list <- list()
for (ii in seq_along(M1M2_scores)) {
  # Modify the title by replacing "_" with a space
  plot_title <- gsub("_", " ", M1M2_scores[ii])
  # Define title color conditionally
  title_color <- if (ii == 1) "red" else "blue"
  # Create the plot
  plot.list[[ii]] <- FeaturePlot(TAM_M1M2_scores, features = M1M2_scores[ii], 
                                 min.cutoff = "q01",
                                 cols = reds) +
    labs(title = plot_title) +
    NoLegend() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 13, face = "plain",
                                colour = title_color)
    )
}
plot_M1M2scores <- ggpubr::ggarrange(plotlist = plot.list, ncol=2, nrow = 1,
                                   common.legend = TRUE, legend = "right")
#
# ============================================================================
# FIGURE 4e: INDIVIDUAL MARKERS
# ============================================================================
Markers <- c("IL1B","CXCL10","TNF","CD80","CD86","IRF1",
             "CD163","MARCO","FN1","IL10","TGFB1","CD274")

# Adding space to labels on edge
new_legends <- paste0("TAM_", single.modules," ") 
names(new_legends) <- single.modules

color_vector <- c(
  "#F8766D",  # A reddish color for FCN1
  "#D89000",  # An orange color for ISG15
  "#A3A500",  # A yellow-green color for CXCL9
  "#39B600",  # A green color for FOSB
  #"#00BF7D",  # A teal color for HLA_II
  "#00BFC4",  # A light blue color for IL1B
  "#619CFF",  # A blue color for FOLR2
  "#F564E3"   # A pink/purple color for APOE
)


# Display median point 
median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}
violin_M1M2 <- list()
for (ii in 1:length(Markers)) {
  #Definiendo color de marcador
  if(Markers[ii]%in%M1_M2$M1_Curated){
    Farbe<-"red"
  }else{
    Farbe<-"blue"
  }
  #Puntos si los genes son poco abundantes
  if (Markers[ii]%in%c("CD80","MARCO","CD274")){
    punto<-0.1
  } else{
    # Remove points for abundant subtypes
    punto<-0
  }
  
  violin_M1M2[[ii]]<-VlnPlot(Poles, features = Markers[ii], pt.size=punto, 
                               assay = "RNA",slot = "data")+
    ylab(Markers[ii])+
    
    #Cambiando orden de x
    scale_x_discrete(limits=single.modules)+
    scale_fill_manual(values = color_vector,
                      labels = new_legends) +
    theme(
      axis.title.y = element_text(size = 10,angle = 90, 
                                  colour = Farbe, vjust = 0.5), 
      #El texto horizontal
      axis.title.x = element_blank(),
      #La etiqueta de cada violin
      axis.text.x = element_blank(),
      #Los numeros de expresion
      axis.text.y = element_text(size = 8),
      #Combinar con otro plot para que no sea grande
      plot.margin=unit(c(-0.5,1,-0.5,1), "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      #El titulo
      plot.title = element_blank(),
      plot.subtitle = element_blank())+
    # Add point at median
    stat_summary(fun = median.stat, geom='point', size = 1, colour = "black")
}
#First part
library(ggpubr)
pre_stack <- ggpubr::ggarrange(plotlist = violin_M1M2, 
                               nrow=length(Markers), ncol = 1,
                             common.legend = TRUE, legend = "right")

Stack_Vln <- annotate_figure(pre_stack, 
                             top = text_grob("M1 and M2 markers",
                                             size = 13, face = "plain"),
                             left = text_grob("Normalized expression", face = "plain", 
                                              color = "black",size = 10, 
                                              hjust =0.5, vjust = 3, rot=90))
#
# ============================================================================
# JOINING
# ============================================================================
C_D <- cowplot::plot_grid(plot_Inf_Repair,plot_M1M2scores,
                          nrow = 2,rel_heights = c(1,0.6))

A_C_D <- cowplot::plot_grid(Plot_MonRes,C_D, ncol = 2)

ABCD <- cowplot::plot_grid(A_C_D, TissMono_Plot, 
                   nrow = 2, rel_heights = c(1,0.7))

cowplot::plot_grid(ABCD,Stack_Vln,
                   ncol = 2, rel_widths = c(1,0.5))
#