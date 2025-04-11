# ============================================================================
# Name: Cluster refinement
# Author: elopez
# Date: May-15-2024
# Description: 
# TODO: 
# ============================================================================
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
# ============================================================================
# Section 1: Trajectory
# ============================================================================
Trajectory <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/MainPhenotype.rds")

#High_cells <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/TissueScores_HighCells_List.rds")

# Extract UMAP dimensions. Select all PCs.
dimred <- Trajectory@reductions$umap@cell.embeddings

# Defining seurat clusters as groups that will be conected (corrected)
# Defining peripheral FCN1 as root of the trajectory
# Root = Blood & FCN1
clustering <- ifelse(Trajectory$MainPheno=="FCN1"&
                     Trajectory$Tissue=="Blood"&
                     Trajectory$seurat_clusters=="C2",
                   "root",
                   as.character(Trajectory$seurat_clusters))

# C3 must be split in CXCL9 and ISG15
clustering <- ifelse(Trajectory$seurat_clusters=="C3"&
                       Trajectory$MainPheno=="CXCL9",
                     paste(Trajectory$seurat_clusters,Trajectory$MainPheno, sep = "_"),
                     clustering)

# Merge C1, C15, and C12 into a single APOE cluster
clustering <- ifelse(clustering%in%c("C1","C12","C15"),"APOE",clustering)

  
library(slingshot)
set.seed(1)
#Proposing blood FCN1 as origin
sds <- slingshot(dimred,
                # Groups
                clustering,
                # Root proposal
                start.clus = "root")

# Obtain linages
sLin <- slingLineages(sds)
sLin
# Save trajectory results
# saveRDS(sds, file = "D:/R_ScriptsPaper/Def_Objects/TissuesSling.rds")
#
# ============================================================================
# Section 2: Plotting
# ============================================================================
#Linages as DF
mst <- slingMST(sds, as.df = TRUE)

#Aesthetic
#Specifying root Monocytes_FCN1
Visualization <- ifelse(Trajectory$MainPheno=="FCN1"&
                          Trajectory$Tissue=="Blood",
                        "Monocytes FCN1",
                        as.character(Trajectory$MainPheno))

# Define the substitutions as a named character vector
#substitutions <- c("FCN1:ISG15" = "CD16 monocytes",
 #                  "HSPs:C1QA" = "Cycling TAMs")
# Perform all the substitutions in a single line using str_replace_all()
#Visualization <- stringr::str_replace_all(Visualization, substitutions)


#Add labels for visualization
df <-data.frame(dimred, "Population" = Visualization)

#Dicrete palette colors
#custom_palette<-DiscretePalette(length(unique(df$Population)),shuffle = T)

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

p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(fill = Population), col = "black", shape = 21, size = 1) +
  # Set the custom combined palette
  scale_fill_manual(values = custom_palette) +
  # Increase the size of the legend circles
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0.5, 
      override.aes = list(size = 3), 
      nrow = 3)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 7)
  )

#Anadir lineas
library(ggrepel)
# Extract the minimum and maximum order for each lineage. 
# In order to know the start and end in each lineage
min_max_data <- mst %>%
  group_by(Lineage) %>%
  summarize(min_order = min(Order), max_order = max(Order))

# Merge the minimum and maximum order information back to the original data
mst <- merge(mst, min_max_data, by = "Lineage")

# Identify whether a point is the "Start" or "End"
mst$point_type <- ifelse(mst$Order == mst$min_order, "Start", ifelse(mst$Order == mst$max_order, "End", "Intermediate"))

# Add red points for "Start" and green points for "End"
p + 
  geom_point(data = mst, aes(color = point_type), size = 3) +
  geom_path(data = mst %>% arrange(Order), aes(group = Lineage), linewidth = 1) +
  scale_color_manual(values = c("Start" = "blue", "End" = "red", "Intermediate" = "black")) +
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
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.title = element_text(size = 7, face = "bold"),
    legend.text = element_text(size = 6)
  )
#