# =============================================================================
# Name: Supplementary Figure S3
# Author: Eric Lopez Huerta
# Date: Jul-21-2024
# Description:
# TODO: Only ggplot2 version 3.4.4 is working with Seurat 
# =============================================================================
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
#CHOOSE DATASET
#Wu, 2021 Dataset
WuTAM_Pheno <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_Pheno.rds")

scores <-colnames(WuTAM_Pheno@meta.data)[grep("score",colnames(WuTAM_Pheno@meta.data))]
Sig.Name <-stringr::str_remove_all(scores, ".score")
# Selecting only the most polarized states
Polarized_TAMs <-subset(WuTAM_Pheno, idents=Sig.Name)

#Azizi, 2018
AziziTAM_Pheno<-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/AziziTAM_Pheno.rds")
# =============================================================================
# SUPP.FIG. S3A: MODULE CO-EXISTANCE BY PAIR FREQUENCY
# =============================================================================
#Frequency of module co-existance
createDotFreqPlot <- function(Modules, AllCombinations, dataset) {
  Frequency <- vector()
  Pairs <- vector()
  for (i in 1:ncol(AllCombinations)) {
    # Regular expression to select strings that have both signatures
    pair.regex <- paste("(?=.*", AllCombinations[1, i], ")", "(?=.*", AllCombinations[2, i], ")", sep = "")
    # Frequency of cells coexpressing the respective pair
    Frequency[i] <- sum(grepl(pair.regex, Modules$Modules, perl = TRUE))
    # Combination name
    Pairs[i] <- paste(AllCombinations[1:2, i], collapse = " & ")
  }
  #Result as DF for plotting
  CombDF <- data.frame(Pairs, Frequency)
  CombDF$Pairs<-stringr::str_replace_all(CombDF$Pairs,"_","-")
  # Order by frequency
  CombDF <- CombDF[order(CombDF$Frequency), ]
  # Fix order for plotting by setting factor and levels
  CombDF$Pairs <- factor(CombDF$Pairs, levels = unique(CombDF$Pairs))
  # Title
  title <- paste0("Module co-existance by pairs ", dataset)
  DotFreq <- ggplot(CombDF, aes(x = Frequency, y = Pairs)) +
    geom_point(colour = "red", size = 2) +
    expand_limits(x = max(Frequency) + 35) +
    ylab("Co-expressing signatures") +
    xlab("Number of cells") +
    geom_text(label = CombDF$Frequency, hjust = -1, size = 2) +
    theme_classic() +
    theme(axis.title = element_text(size = 8),
          axis.text.y = element_text(size = 5),
          axis.text.x = element_text(size = 6),
          # Adjust size and center
          plot.title = element_blank(),
          panel.border = element_rect(fill = NA)) +
    ggtitle(title)
  
  return(DotFreq)
}

# Create the AllCombinations matrix
All.Combinations <- combn(unique(Idents(Polarized_TAMs)), 2)

# Apply the function to both WuTAM_Pheno and Azizi_Modules and store the results in a list
dotfreq_plots <- list(
  Wu_plot = createDotFreqPlot(WuTAM_Pheno, All.Combinations, "(Wu, 2021)"),
  Azizi_plot = createDotFreqPlot(AziziTAM_Pheno, All.Combinations, "(Azizi, 2018)")
)

# Now, to print the plots, you can use
#Two.Dots.Plots <- cowplot::plot_grid(plotlist= dotfreq_plots, nrow = 2)
# =============================================================================
# SUPP.FIG. S3B: HEATMAP HCLUST
# =============================================================================
#Wu, 2021 enrichment DF
Wu_enrich <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_ES_DF_TAM.rds")
Wu_enrich <- Wu_enrich[,c(1:3,6,7)]

#Azizi, 2018 enrichment DF
Azizi_enrich <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Azizi_ES_DF_TAM.rds")

library(pheatmap)   
library(gplots)
#library(ComplexHeatmap)
#library(grid)

HeatClust <- function(DF_enrich, dataset) {
  # Removing "_"
  colnames(DF_enrich) <- stringr::str_replace_all(colnames(DF_enrich), "_", " ")
  rownames(DF_enrich) <- NULL
  
  # Palette color
  Morado <- RColorBrewer::brewer.pal(9, "Purples")
  Morado <- Morado[2:9]
  
  # Title
  paste_title <- paste0(nrow(DF_enrich), " TAMs ", dataset)
  
  # Suppress automatic plotting with print = FALSE
  C <- pheatmap(t(DF_enrich), col = Morado, cluster_cols = TRUE, cluster_rows = TRUE,
                #fontsize_col = 0, fontsize_row = 6,
                fontsize = 6,
                border_color = NA,
                #main = paste_title,
                annotation_legend = FALSE, annotation_names_row = FALSE
                )
  
  # Convert to ggplot for use with patchwork
  Heat_morado <- ggplotify::as.ggplot(C)
  return(Heat_morado)
}

# Apply the function to both Wu_enrich and Azizi_enrich and store the results in a list
heatmap_scores <- list(
  Wu_plot = HeatClust(Wu_enrich, "(Wu, 2021)"),
  Azizi_plot = HeatClust(Azizi_enrich, "(Azizi, 2018)")
)
#corroborate plot
#Two.Heat.Plots <- cowplot::plot_grid(plotlist= heatmap_scores, nrow = 2)
# =============================================================================
# SUPP.FIG. S3C: AVERAGE ENRICHMENT LONG HEATMAP 
# =============================================================================
#Select columns with score and phenotype data
ES_Phen <- grep("score|Modules", colnames(WuTAM_Pheno@meta.data))
#Make a df for each phenotype
DFL<-dplyr::group_split(WuTAM_Pheno@meta.data[,ES_Phen], Modules)

#Calculate average ES for each phenotype
avg.ES<-function(df){
  #Excluding 7th column
  means<-colMeans(df[,-ncol(df),drop=F])
  #Number of cells for each phenotype
  num_rows<-nrow(df)
  return(c(means,Labeled_cells=num_rows))
}
meanES_list<-lapply(DFL, avg.ES)
#Storing name of the respective phenotype in vector
phen.nam<-vector()
for (i in 1:length(DFL)) {
  phen.nam[i]<-unique(DFL[[i]]$Modules)
}
#Name each element of the list
names(meanES_list)<-phen.nam

#Convert list into DF
proms<-plyr::ldply(meanES_list,rbind)
#Adding number of cells to the phenotype name to store in Excel
proms$PhenNumCells<-paste(proms$.id," (n=",proms$Labeled_cells,")",sep = "")

#Order by abundance, but leaving "Unassigned" at the bottom.
prom.sorted<-proms[!proms$.id=="Unassigned",]
prom.sorted<-prom.sorted[order(prom.sorted$Labeled_cells, decreasing = T),]
#Joining bottom containing Unassigned
prom.sorted<-rbind(prom.sorted,proms[proms$.id=="Unassigned",])
#write.csv(prom.sorted, file = "D:/R_ScriptsPaper/Def_Objects/PROMEDIO_ES.csv", row.names = F)

#Only columns that will be plotted in long format
DF_COMPLETE<-reshape2::melt(prom.sorted[,!colnames(prom.sorted)%in%c("Labeled_cells","PhenNumCells")])  
#Renaming columns to have consistency
colnames(DF_COMPLETE)<-c("Label","score","value")
#Aesthetic: Removing "_"
DF_COMPLETE$score<-stringr::str_replace_all(DF_COMPLETE$score,"_score"," score")

#VISUALIZE TOP DOBLE COMBINATIONS
print(head(prom.sorted$.id[stringr::str_count(prom.sorted$.id,":")==1]),10)
#CHOSEN BASED ON FUNCTIONS:
Dobles<-c("ISG15:CXCL9","FOLR2:APOE","FCN1:IL1B")

#VISUALIZE TOP TRIPLE COMBINATIONS
print(head(prom.sorted$.id[stringr::str_count(prom.sorted$.id,":")==2]),10)
#CHOSEN: 
#Triples<-c("FOSB:FOLR2:APOE","FCN1:ISG15:IL1B")

#ORDEN ESCALONADO
# Define the substitutions as a named character vector
# Names are the second values
substitutions <- c("FCN1" = "a_FCN1",
                   "IL1B" = "b_IL1B",
                   "ISG15" = "c_ISG15",
                   "CXCL9" = "d_CXCL9",
                   "FOSB" = "e_FOSB",
                   "FOLR2" = "f_FOLR2",
                   "APOE" = "g_APOE")

#Working only with chosen phenotypes
chosen<-c(names(substitutions),Dobles,
          #Triples,
          "Unassigned")
dfggplot<-DF_COMPLETE[DF_COMPLETE$Label%in%chosen,]
# Perform all the substitutions to set order
a_substituted <- unique(stringr::str_replace_all(dfggplot$Label, substitutions))
#Ordenado con prefijo alfabetico
a_subsort<-rev(sort(a_substituted))
cleaned_strings <- gsub("[a-h]_", "", a_subsort)
Yorder <- cleaned_strings
#Order in x axis
Xorder <- paste0(names(substitutions)," score")
#Paleta de colores azul para -1 y rojo para 1
paletaAzul<-rev(RColorBrewer::brewer.pal(n = 11, name = 'RdBu'))


LongESHeat <- ggplot(dfggplot, aes(x = score, y = Label, fill = value)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradientn(colours = paletaAzul, limits = c(-1, 1)) +
  scale_y_discrete(limits = Yorder) +
  scale_x_discrete(limits = Xorder, position = "top") +
  ylab("TAM phenotype") + xlab(NULL) + 
  labs(fill = "") +
  geom_text(label = round(dfggplot$value, 2), size = 3.5) +
  theme(
    axis.text.x = element_text(angle = 320, colour = "red", size = 11),
    axis.text.y = element_text(size = 11),
    axis.text.x.top = element_text(vjust = 0.1, hjust = 1),
    axis.title.y = element_text(size = 13),
    legend.title = element_blank(),
    #legend.key.size = unit(0.4, "cm"),  # Adjust the key size
    #legend.key.height = unit(0.3, "cm"),  # Adjust the key height
    #legend.key.width = unit(0.5, "cm"),  # Adjust the key width
    # Adjust the margin
    #legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),  
    #legend.spacing = unit(0.2, "cm"),  # Adjust the spacing between legend items
    legend.position = "top",  # Position legend at the top
    legend.justification = "left",  # Center the legend horizontally
    legend.direction = "horizontal",  # Arrange legend items horizontally
    panel.border = element_rect(fill = NA),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank()
  )
# =============================================================================
# SUPP.FIG. S3E: PATHWAYS LONG
# =============================================================================
GSVA_Res <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_gsvaPaths.rds")

#Adding enrichment data to seurat object
TAM_paths <-AddMetaData(WuTAM_Pheno, metada=GSVA_Res)

# In order to calculate average ES, we will split df by identity/Modules column
Idents(TAM_paths) <- TAM_paths$Modules

# Average ES in polarized pops
S.Obj <-subset(TAM_paths, idents =c(Sig.Name,Dobles))

# Only columns that will be used
ES_DF <-  S.Obj@meta.data[,grep("KEGG|GOBP|HALLMARK|Modules",colnames(S.Obj@meta.data))]
DFL <- dplyr::group_split(ES_DF, Modules)
mean_ES <- list()
mod_ES <- vector()

for (gg in 1:length(DFL)){
  #Selecting only numeric values
  d.f <-as.data.frame(DFL[[gg]][-1])
  mean_ES[[gg]] <-colSums(d.f)/nrow(d.f)
  mod_ES[gg] <-as.character(unique(DFL[[gg]]$Modules))
}

names(mean_ES) <- mod_ES
# List to DF
promsPaths<-plyr::ldply(mean_ES,rbind)

#---------------------GROUPING BY REDUNDANCY
library(stringdist)
library(dplyr)
library(stringr)

# Function to remove redundant terms and include the highest value for each pathway
remove_redundant_terms <- function(terms, avg_ES, threshold, n.rows) {
  to_eliminate <- "GOBP_|KEGG_|HALLMARK_|POSITIVE_|NEGATIVE_|REGULATION_OF_|DEFENSE_|RESPONSE_TO_|_RESPONSE$|SIGNALING_PATHWAY|_PROCESS$|CELLULAR_|_CELL$|CELL_|^TRANS"
  # Clean up the terms by removing redundant phrases
  cleaned_terms <- str_remove_all(terms, to_eliminate)
  cleaned_terms <- str_replace_all(cleaned_terms, "_", " ")
  
  # Compute a distance matrix based on string similarity
  # Jaro-Winkler similarity
  dist_matrix <- stringdistmatrix(cleaned_terms, cleaned_terms, method = "jw")  
  
  # Group terms based on similarity threshold
  # Higher threshold for more flexible clustering
  groups <- hclust(as.dist(dist_matrix)) %>%
    cutree(h = threshold)
  
  # Find the highest value for each pathway
  # Get the max value for each column
  highest_values <- apply(avg_ES[,-1], 2, max, na.rm = TRUE)
  
  # Combine terms and their groups into a dataframe
  df_1 <- data.frame(
    row.names = NULL,
    Term = terms,
    Group = groups,
    HighestValue = highest_values
  )
  
  # Arrange by the highest value
  # Select only the first row of each group with the smallest P.DE
  df_2 <- df_1 %>%
    arrange(Group, desc(HighestValue)) %>%  # Sort by group and descending HighestValue
    group_by(Group) %>%
    # Conditionally select rows: 2 if group size >= 3, 1 otherwise
    filter(n() >= 3 & row_number() <= n.rows | n() < 3 & row_number() == 1) %>%
    ungroup()
  
  return(df_2)
}

# Pathway names. Exclude the first column (identifiers)
terms <- colnames(promsPaths)[-1]
reduced_terms <- remove_redundant_terms(terms=terms, avg_ES=promsPaths, 
                                        threshold = .3, n.rows = 2)

# Display only reduced terms
short_proms<-promsPaths[,c(".id", reduced_terms$Term)]

# using abbreviations
abreviations <- c("ANTIGEN_PROCESSING_AND" = "AG",
                  "NEGATIVE_REGULATION_OF" = "(-)",
                  "POSITIVE_REGULATION_OF" = "(+)",
                  "REGULATION_OF" = "REG",
                  #"CELLULAR_" = "",
                  "IMMUNOGLOBULIN" = "IG",
                  "HALLMARK" = "MSigDB")
#
abrev_paths <- stringr::str_replace_all(colnames(short_proms), abreviations)
colnames(short_proms)<-stringr::str_replace_all(abrev_paths,"_"," ")

#pheatmap
dat <- tibble::column_to_rownames(short_proms, var = ".id")

library(ComplexHeatmap)
library(grid)
reds.pal <- RColorBrewer::brewer.pal(9, "YlOrRd")
# Create the heatmap
ht <- Heatmap(t(dat), 
              col = reds.pal, 
              name = "NES",
              cluster_columns = TRUE, 
              cluster_rows = TRUE, 
              show_column_names = TRUE, 
              show_row_names = TRUE,
              column_names_gp = gpar(fontsize = 7.5), 
              row_names_gp = gpar(fontsize = 5),
              width = unit(5, "cm"),  # Adjust the heatmap width
              heatmap_legend_param = list(
                title = "",
                title_gp = gpar(fontsize = 8, fontface = "bold"),
                labels_gp = gpar(fontsize = 6),
                legend_height = unit(3, "cm"),
                #legend_position = "bottom",
                legend_direction = "horizontal"
                ))

# Convert the ComplexHeatmap to a ggplot object and plot legend on top
LongPaths <- ggplotify::as.ggplot(~draw(ht, heatmap_legend_side = "top"))
# =============================================================================
# DISPLAYING PLOTS
# =============================================================================
# Upper part
jpeg("D:/R_ScriptsPaper/Figures/Manuscrito/pre_SuppFigS3A-B.png", 
     res = 600,width = 18, height = 5, units = "cm")
cowplot::plot_grid(dotfreq_plots[[1]],heatmap_scores[[1]],
                   ncol = 2)
dev.off()

# Lower part
cowplot::plot_grid(LongESHeat,NULL,LongPaths,ncol = 3,
                   rel_widths = c(0.8,0.6,1.2))
