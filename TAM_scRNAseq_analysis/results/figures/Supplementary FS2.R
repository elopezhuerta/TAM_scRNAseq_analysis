# =============================================================================
# Name: Supplementary Figure S2
# Author: Eric Lopez Huerta
# Date: Jul-14-2024
# Description: 
# TODO: Import ggplot2 version 3.4.4 before Seurat to avoid error 
# =============================================================================
library(ggplot2)
library(Seurat)
library(ggpubr)
library(dplyr)
library(gridExtra)
#Wu DATASET
TAMs <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAMs.rds")

# Setting clusters as idents
Idents(TAMs) <-paste0("C",TAMs$seurat_clusters)

#Wu DATASET
Wu_SortGenes <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_SortGenes.rds")
Modules <- Wu_SortGenes
#
# =============================================================================
# SUPP.FIG. S2A: OVERREPRESENTED PATHS
# =============================================================================
# Import over-representation results 
go_kegg_DFresults <- read.csv(file = "D:/R_ScriptsPaper/Def_Objects/Wu_go_kegg_DFresults.csv", 
                              header = T)

# Setting order
go_kegg_DFresults$Module <- factor(go_kegg_DFresults$Module ,
                                   levels = c("FCN1","ISG15","CXCL9","IL1B","FOSB","FOLR2","APOE","MKI67","CALR","ACTB"))

# Filter and choose number of pathways to represent 
go_kegg_filtered <- go_kegg_DFresults%>%
  # Group by Ont and globally arrange by P.DE and DE within Ont
  group_by(Ont, Module) %>%
  arrange(P.DE, desc(DE), .by_group = TRUE) %>%
  # Extract only the top pathways per subtype
  slice_head(n= 6) %>%
  ungroup()

# Shorten long names
go_kegg_filtered %>%
  mutate(LetterSpaceCount = nchar(gsub("[^a-zA-Z ]", "", Term))) %>%  # Count letters and spaces
  filter(LetterSpaceCount >= 60)%>%
  pull(Term)

library(stringr)
go_kegg_filtered$Term<-str_replace_all(go_kegg_filtered$Term,
                                       "biological process involved in | cytokine and ",
                                       " ")
# Separating GO results from KEGG
Ont_df_list <- group_split(go_kegg_filtered, Ont)

# Plotting
library(ggplot2)
library(cowplot)

# Colour palette
greens.pal <- RColorBrewer::brewer.pal(9, "Greens")[2:9]

# Generate plot for each database
base <- c("GO", "KEGG")
PLIST <- list()
for (ii in seq_along(Ont_df_list)) {
  DF <- as.data.frame(Ont_df_list[[ii]])
  DF$Pathway <- factor(DF$Term, levels = unique(DF$Term))
  DF$Module <- factor(DF$Module, levels = unique(DF$Module))
  DF$GeneRatio <- DF$DE / DF$N
  
  # Define the positions for vertical dashed lines
  midpoints <- seq(0.5, nrow(DF) + 0.5, by = 1)
  
  PLIST[[ii]]<-ggplot(DF, aes(x= Module, y= Pathway, fill=GeneRatio))+
    geom_tile()+
    #setting white background
    theme_classic()+
    #Choosing colour and title
    scale_fill_gradientn(colours = greens.pal)+
    ggtitle(paste(base[ii]," pathways"))+
    #x thick marks orientation and vertical dashlines 
    theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.8)),
          panel.border = element_rect(fill=NA),
          axis.title = element_blank(),
          axis.text.y = element_text(size = rel(0.8)),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          # Adjust margins to avoid truncation
          plot.margin = margin(10, 10, 10, 20)  # Top, Right, Bottom, Left
    ) +
    # Add vertical dashed lines
    geom_vline(xintercept = midpoints, linetype = "dashed")
  
}
# Plot side by side
Plot_go_kegg <- cowplot::plot_grid(plotlist=PLIST, 
                                   ncol = 2)
# =============================================================================
# SUPP.FIG. S2B: INTERSECTING GENES BETWEEN WU AND AZIZI MODULES
# =============================================================================
Wu_fragmented <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_fragmented.rds")
# Selecting and renaming
Wu_fragmented <- Wu_fragmented[c(2,1,4,8,9,6,7)]
names(Wu_fragmented) <-paste0("Wu_", c("FCN1","ISG15","CXCL9","IL1B","FOSB","FOLR2","APOE"))

Azizi_fragmented <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Azizi_fragmented.rds")
# Adding prefix to avoid confusion
names(Azizi_fragmented) <-paste0("Azizi_", c("FOLR2","APOE","HLA_II","FCN1","ISG15"))

# Function to count matches
count_matches <- function(x, y) {
  sum(sapply(y, function(pattern) sum(x %in% pattern)))
}

# Apply the function to each element of Azizi_fragmented using lapply
match_counts <- lapply(Azizi_fragmented, function(sub_list) {
  sapply(Wu_fragmented, count_matches, y = sub_list)
})

# Convert the match_counts to a data frame
df <- as.data.frame(match_counts)

# Reshape the data frame to long format
df <-tibble::rownames_to_column(df,var = "x_axis")

# Reshape the data frame to long format
df_long <- tidyr::pivot_longer(df, cols = -x_axis, 
                               names_to = "y_axis", values_to = "Coincidences")

# Calculate the percentage of coincidences based on the size of the corresponding module in Azizi_fragmented
df_long <- df_long %>%
  mutate(Percentage = ifelse(y_axis %in% names(Azizi_fragmented), 
                             Coincidences / sapply(y_axis, function(x) length(Azizi_fragmented[[x]])) * 100,
                             NA))


# Setting X order by percentage
Xorder <-unique(df_long$x_axis[order(df_long$Percentage, decreasing = T)])

df_long$x_axis <- factor(df_long$x_axis, 
                         levels = Xorder)

# Setting Y order by percentage
Yorder <-unique(df_long$y_axis[order(df_long$Percentage, decreasing = T)])
# aesthetic order
Yorder <- Yorder[c(1:3,5,4)]
df_long$y_axis <- factor(df_long$y_axis,
                         levels= rev(Yorder))

#HEATMAP
spect.pal<-RColorBrewer::brewer.pal(9,"Spectral")

Plot_Coincidences <- ggplot(df_long, 
                            aes(x = x_axis, y = y_axis, fill = Percentage)) +
  geom_tile() +
  theme_classic() +
  # Colour palette
  scale_fill_gradientn(colours = rev(spect.pal)) +
  # Label order by percentage
  scale_x_discrete(position = "top") +
  # Axis titles
  ylab("Dataset: Azizi,2018") + xlab("Dataset: Wu,2021") +
  labs(fill = "% of Azizi modules") +
  # Number inside the heatmap squares with % format
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), size = 2.5) + 
  
  theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5,
                                   size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 8),
        panel.border = element_rect(fill = NA),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank())

#
# =============================================================================
# ARRENGING ALL PLOTS
# =============================================================================
Plot_go_kegg

# Second part
jpeg("D:/R_ScriptsPaper/Figures/Manuscrito/pre_SuppFigS2B.png", 
     res = 600,width = 18, height = 5.5, units = "cm")

Plot_Coincidences

dev.off()
#