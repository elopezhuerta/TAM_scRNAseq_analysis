# =============================================================================
# Name: Supplementary Figure S1
# Author: Eric Lopez Huerta
# Date: Jul-14-2024
# Description: 
# TODO:
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
# SUPP.FIG. S1B-C: DISTRIBUTION OF MODULE GENES ACROSS CLUSTERS
# =============================================================================
# Analysis on large modules
large_modules <- Modules[lengths(Modules) > 100]

# Marker of each cluster
TAMDF <-FindAllMarkers(TAMs,features = as.character(unlist(large_modules)),
                       # Calculate on normalize data
                       assay = "RNA", slot = "data", 
                       # Minimal fold change and expression %
                       logfc.threshold = 0.4, 
                       min.pct = 0.33, 
                       only.pos=T)

#-------------------COINCIDENCES
top_gene_filtered <- TAMDF%>%
  filter(avg_log2FC >= 0.4, p_val_adj<0.01)
# Get as list
DEGs_clusters <- split(top_gene_filtered$gene, top_gene_filtered$cluster)

# Function to count matches between DEGs and large modules
count_matches <- function(x, y) {
  sum(sapply(y, function(pattern) sum(x %in% pattern)))
}

# Apply the function to each element of DEGs_clusters using lapply
match_counts <- lapply(DEGs_clusters, function(sub_list) {
  sapply(large_modules, count_matches, y = sub_list)
})
# Calculate non-DEGs for each module
non_DEGs_counts <- lengths(lapply(large_modules, setdiff, unlist(DEGs_clusters)))

# Convert the match_counts to a data frame
df_DEGs <- as.data.frame(match_counts)
# Reshape the data frame to long format
df_DEGs <- tibble::rownames_to_column(df_DEGs,var = "Modules")

# Convert non-DEGs to a data frame
df_non_DEGs <- data.frame(
  Modules = names(non_DEGs_counts),
  non_DEGs = unname(non_DEGs_counts)
)

# Join the data frames
df_combined <- df_DEGs %>%
  left_join(df_non_DEGs, by = "Modules")

df_long <- tidyr::pivot_longer(df_combined, cols = -Modules, 
                               names_to = "DEGs", values_to = "Coincidences")

# Calculating % of coincides based on modules
df_long <- df_long %>%
  mutate(Percentage = ifelse(Modules %in% names(large_modules), 
                             Coincidences / sapply(Modules, function(x) length(large_modules[[x]])) * 100,
                             NA)) %>%
  # Setting X and Y
  arrange(desc(Percentage))


module_plots <- list()
# This will also dictate order
my_module <- c("turquoise", "blue")
for (module in my_module) {
  # Filter data for the current module
  module_data <- df_long %>%
    filter(Modules == module) %>%
    # This ensures that non_DEGs column is placed at the end
    mutate(DEGs = factor(DEGs, levels = c(setdiff(DEGs, "non_DEGs"), "non_DEGs")))
  
  # Create the plot for the current module
  plot <- ggplot(module_data, aes(x = DEGs, y = Percentage, fill = DEGs)) +
    geom_bar(stat = "identity", color = "black") +
    geom_text(aes(label = Coincidences), vjust = -0.5, size = 3) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(0, 80),
      breaks = seq(0, 80, by = 10)
    ) +
    scale_fill_manual(
      values = c(rep(module, length(levels(module_data$DEGs)) - 1), "gray")
    )+
    labs(
      title = paste(module, "module","(n=",length(large_modules[[module]]),"genes)"),
      x = "DEGs per Seurat cluster",
      y = "Module Percentage"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  module_plots[[module]] <- plot
}

Plot_DEGs_Clusters <- patchwork::wrap_plots(module_plots, nrow = 1, ncol = 2)
# Proportion
cowplot::plot_grid(Plot_DEGs_Clusters, NULL, nrow = 2, rel_heights = c(1,0.5))
#
#