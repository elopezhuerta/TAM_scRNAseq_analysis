# ============================================================================
# Name: Figure 1-PhenoPaths
# Author: elopez
# Date: Jun-30-2024
# Description: 
# TODO: 
# ============================================================================
library(Seurat)
library(ggplot2)
library(ggpubr)
WuTAM_Pheno <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_Pheno.rds")
# Add C prefix
WuTAM_Pheno$seurat_clusters <-paste0("C",WuTAM_Pheno$seurat_clusters)
# Colores de UMAP
reds.pal <-RColorBrewer::brewer.pal(9, "YlOrRd")
spectral.pal <-rev(RColorBrewer::brewer.pal(11,"Spectral"))
# Extracting scores
scores <-colnames(WuTAM_Pheno@meta.data)[grepl(".score",colnames(WuTAM_Pheno@meta.data))]
#Polarized states
Polarized <- stringr::str_remove_all(scores,".score")
#TAM labels (altered)
TAM_Labs<-paste0("TAM_",Polarized)
#
# ============================================================================
# Figure 1A: PHENOTYPE PER CLUSTER
# ============================================================================
#ASSIGNING MOST ABUNDANT PHENOTYPE PER CLUSTER TO VISUALIZE
all_cl <- table(WuTAM_Pheno$seurat_clusters)

#FOR SECOND ROUND
#List of cells with high ES
HighCells_List <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuScores_HighCells_List.rds")
#HighCells_List$Unassigned <- setdiff(colnames(WuTAM_Pheno),unique(unlist(HighCells_List)))

#Transform list into named vector
my_vector <-unlist(HighCells_List)
names(my_vector) <-rep(names(HighCells_List), lengths(HighCells_List))

# Initialize placeholders
substitutions <- vector()
overlap.map <- list()
# Define threshold to name cluster
threshold <- 33.33
for (ii in seq_along(all_cl)) {
  # Get the names that are in the specified cluster
  matched_names <- names(my_vector)[
    my_vector %in% colnames(WuTAM_Pheno)[WuTAM_Pheno$seurat_clusters == names(all_cl[ii])]]
  # Calculate % that each phenotype occupies in each cluster
  percent <- table(matched_names) / all_cl[ii] * 100
  sorted_values <- percent[order(percent, decreasing = TRUE)]
  sorted_names <- names(sorted_values)
  # The Highest phenotype should be at least 1/3 of the cluster
  if (length(sorted_values) > 0 && sorted_values[1] >= threshold) {
    substitutions[ii] <- sorted_names[1]
    # If the second higher is also above 33.33%, create a list
    # with the respective cluster and the name of the second highest phenotype as element
    if (length(sorted_values) > 1 && sorted_values[2] >= threshold) {
      overlap.map[[names(all_cl[ii])]] <- sorted_names[2]
    }
  } else {
    # If no phenotype achieve parameters. Label as unassigned
    substitutions[ii] <- "Unassigned"
    # Rescue in unassigned
    overlap.map[[names(all_cl[ii])]] <- sorted_names[1]
    # Overlap with the highest of 2 highest if difference is less than 5%
    if (sorted_values[1]-sorted_values[2] <=5 ){
      overlap.map[[names(all_cl[ii])]] <-c(sorted_names[1],sorted_names[2])
    }
  }
}
names(substitutions) <- names(all_cl)

# Replace Seurat clusters with most abundant phenotype and store in metadata
WuTAM_Pheno$MainPheno <- dplyr::recode(WuTAM_Pheno$seurat_clusters, !!!substitutions)

# RESCUE PHENOTYPES
# Rescue phenotypes (as names) and set respective overlaping phenos
for (ii in seq_along(overlap.map)) {
  # Cells within the clusters
  cl_logic<-WuTAM_Pheno$seurat_clusters%in%names(overlap.map[ii])
  # Cells with high of the respective phenotype
  High_logic<-colnames(WuTAM_Pheno)%in%as.character(unlist(HighCells_List[overlap.map[[ii]]]))
  # Cells with both
  logical_vect <-cl_logic&High_logic
  # Substituing
  WuTAM_Pheno$MainPheno<-ifelse(logical_vect,
                           overlap.map[[ii]],
                           as.character(WuTAM_Pheno$MainPheno))
}

#AESTHETIC
WuTAM_Pheno.plot <- WuTAM_Pheno
#Add TAM prefix
WuTAM_Pheno.plot$MainPheno <-paste0("TAM_",WuTAM_Pheno.plot$MainPheno)

#Run until you find a palette you like
#Distinct<-DiscretePalette(length(unique(WuTAM_Pheno.plot$MainPheno)),shuffle = T)
#Chosen color palette 
Distinct <- c("#FFE100", #APOE
              "#0075DC", #CXCL9
              "#2BCE48", #FCN1
              "#FF0010", #FOLR2
              "#993F00", #FOSB
              "#FF5005", #IL1B
              "#F0A0FF", #ISG15
              "darkgrey" #Unassigned
)

ClPerPheno <- DimPlot(WuTAM_Pheno.plot, group.by = "MainPheno", cols = Distinct)+
  labs(title = "TAM subsets")+
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "plain"),
        axis.title = element_text(size = 10),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.spacing.x= unit(0.1, 'mm'),
        legend.spacing.y= unit(0.1, 'mm'))+
  guides(color=guide_legend(nrow=4, byrow=TRUE,
                            override.aes = list(size=3)))
#
# ============================================================================
# Figure 1B: SCORE PLOTS
# ============================================================================
# Rename those columns by removing underscores
colnames(WuTAM_Pheno@meta.data)[colnames(WuTAM_Pheno@meta.data) %in% scores] <- 
  gsub("_", " ", scores)

PhenoScores <- FeaturePlot(WuTAM_Pheno, features = gsub("_", " ", scores), ncol = 4,
                          max.cutoff = "q99",
                          cols = spectral.pal) &
  NoLegend() &
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "plain"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

library(patchwork)
PhenoScores <- PhenoScores + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right")
#
# ============================================================================
# Figure 1C: SIGNATURE WITH DOTPLOTS
# ============================================================================
#Plot only most polzarized states
Poles <- subset(WuTAM_Pheno, idents=Polarized)
#Cluster order
levels(Poles) <- Polarized

#Signature genes to plot
#Intersecting signatures + rescued modules
BestGenes <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_BestGenes.rds")

# Top ten genes
BG <-lapply(BestGenes, head,5)
# Ensuring that best genes are displayed for each signature
BG$CXCL9[5] <- "CXCL11"
BG$IL1B[c(1,3,4)] <- c("IL1B","PNRC1","BTG1")
BG$FOSB [3:5] <- c("HSPA1B", "JUN", "DUSP1")
BG$FOLR2[c(4,5)] <- c("FOLR2", "MRC1")
BG$APOE <- c("APOE","APOC1", "CTSD", "CTSB", "PSAP")

# Dotplot
SignPattern_Dots <- DotPlot(WuTAM_Pheno, features = as.character(unlist(BG)), 
                            assay = "RNA", cols = "Spectral", 
                            dot.scale = 3) + # Reduce dot size with dot.scale
  scale_y_discrete(limits = 
                     #Setting order
                     rev(Polarized),
                   labels = 
                     #Changing labels
                     rev(TAM_Labs)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 10),
        axis.title = element_blank(),
        legend.title = element_text(size = 9),  # Adjust legend title size
        legend.text = element_text(size = 9))   # Adjust legend text size
#
# ============================================================================
# Figure 1D: PATHWAYS SUMMARY HEATMAP
# ============================================================================
Poles <- subset(WuTAM_Pheno, idents = Polarized)

# List of pathways of interest
CustomPaths <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/CustomPaths.rds")

GSVA_Res <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_gsvaPaths.rds")

# Adding enrichment data to seurat object
TAM_Custom <- AddMetaData(Poles, metada=GSVA_Res[,CustomPaths])
Idents(TAM_Custom) <- TAM_Custom$Modules

library(dplyr)
# Only representative Pathways and Modules column
ES_DF <-TAM_Custom@meta.data[,c("Modules",CustomPaths)]

DFL<-dplyr::group_split(ES_DF, Modules)
mean_ES<-list()
mod_ES<-vector()
for (gg in 1:length(DFL)){
  d.f<-as.data.frame(DFL[[gg]][-1])
  mean_ES[[gg]]<-colSums(d.f)/nrow(d.f)
  #Creando vector con el nombre
  mod_ES[gg]<-as.character(unique(DFL[[gg]]$Modules))
}
names(mean_ES) <-mod_ES
# List to DF
promsPaths <-plyr::ldply(mean_ES,rbind)

#pheatmap
library(tibble)
dat <-column_to_rownames(promsPaths, var = ".id")
#Setting order in X (rownames)
dat <-dat[Polarized,]

# Aesthethic
colnames(dat) <-stringr::str_replace_all(colnames(dat),"_"," ")
colnames(dat) <-stringr::str_replace_all(colnames(dat),"HALLMARK","MSIGDB")

# Changing names long names
short_names <- c("GOBP ANTIGEN PRESENTATION",
                 "GOBP POSITIVE REGULATION ERK1/2 CASCADE",
                 "GOBP LIPOPROTEIN REGULATION")
colnames(dat)[grep("ENDOGEN|LEVELS|ERK",colnames(dat))] <- short_names

# Re order for plotting
dat <- dat[,c(1:3,5,8,4,6:7,9:ncol(dat))]

# Add prefix
rownames(dat) <- paste("TAM",rownames(dat),sep = "_")

library(pheatmap)   
library(gplots)

S <- pheatmap(t(dat), col=spectral.pal, cluster_cols=F, cluster_rows=F, 
            fontsize_col=6, fontsize_row = 6, legend_labels = "NES",
            border_color=NA, angle_col = 45)

# Turn into a ggplot to use it with patchwork
ResumPaths<-ggplotify::as.ggplot(S)
#
# ============================================================================
# ARRENGING PLOTS (Figure 1A-D: )
# ============================================================================
AB <- cowplot::plot_grid(ClPerPheno,PhenoScores, ncol = 2, 
                        #labels = c("a","b"), label_size = 19,
                        rel_widths = c(0.5,1))
CD <- cowplot::plot_grid(SignPattern_Dots, ResumPaths, ncol = 2,
                        #labels = c("c"), label_size = 19,
                        rel_widths = c(1,0.6))

# Final
cowplot::plot_grid(AB,CD, nrow = 2,rel_heights = c(1,0.6))
