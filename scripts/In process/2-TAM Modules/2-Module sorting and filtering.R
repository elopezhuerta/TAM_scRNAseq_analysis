# =============================================================================
# Name: Module sorting and filtering
# Author: elopez
# Date: May-03-2024
# Description: Sorting genes by membership and determining if module distribution

# TODO: 
# =============================================================================
#GENERAL USE
library(Seurat)
library(dplyr)
library(ggplot2)
#Palette
reds.pal<-RColorBrewer::brewer.pal(9, "YlOrRd")
greens.pal<-RColorBrewer::brewer.pal(8,'Greens')
spect.pal<-RColorBrewer::brewer.pal(11,"Spectral")

#CHOOSE DATASET TO WORK WITH:

#Wu DATASET
TAMs <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAMs.rds")

# AZIZI DATASET
TAMs <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Azizi_TAMs.rds")

# Setting clusters as idents
Idents(TAMs) <-paste0("C",TAMs$seurat_clusters)
#
# =============================================================================
# Section 1: SORTING AND DEFINING NEGATIVE MODULES
# =============================================================================
# Function to sort modules by membership.
# Placing genes with negative values in another modules with prefix anti_
SortModuleList <-function(module.list, memberships, min.size){
  # Remove grey genes
  module.list["grey"]<-NULL
  # Initialize a list to store sorted genes, including "anti_" modules
  SortGenes <- list()
  # Loop through each module in module.list
  for (module in names(module.list)) {
    # Subset memberships for the current module's genes
    dat_module <- memberships %>%
      filter(rownames(memberships) %in% module.list[[module]])
    
    # Separate positive and negative values for the respective module
    positive_genes <- dat_module %>%
      filter(!!sym(module) >= 0) %>%
      arrange(desc(!!sym(module))) %>%
      rownames()
    
    negative_genes <- dat_module %>%
      filter(!!sym(module) < 0) %>%
      # Sort in ascending order for negatives
      arrange(!!sym(module)) %>%  
      rownames()
    # Store positive genes under the original module
    SortGenes[[module]] <- positive_genes
    # Store negative genes under the "anti_" module
    if (length(negative_genes) > min.size) {
      SortGenes[[paste0("anti_", module)]] <- negative_genes
    }
  }
  
  return(SortGenes)
}

# Wu DATASET
Wu_ModuleGenes<-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_ModuleGenes.rds")
Wu_datKME<-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_datKME.rds")
# Sort Wu modules by membership
Wu_SortGenes <-SortModuleList(module.list = Wu_ModuleGenes,
                              memberships = Wu_datKME,
                              min.size = 30)

# saveRDS(Wu_SortGenes, file = "D:/R_ScriptsPaper/Def_Objects/Wu_SortGenes.rds")


# AZIZI DATASET
Azizi_ModuleGenes <-readRDS("D:/R_ScriptsPaper/Def_Objects/Azizi_ModuleGenes.rds")
Azizi_datKME<-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Azizi_datKME.rds")
# Sort Wu modules by membership
Azizi_SortGenes <-SortModuleList(module.list = Azizi_ModuleGenes,
                                 memberships = Azizi_datKME,
                                 min.size = 15)

# saveRDS(Azizi_SortGenes, file = "D:/R_ScriptsPaper/Def_Objects/Azizi_SortGenes.rds")
# ============================================================================
# Section 2: FRAGMENTING LARGE MODULES WITH DEGs (WU DATASET)
# ============================================================================
#Wu DATASET
Wu_SortGenes <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_SortGenes.rds")
Modules <- Wu_SortGenes

# Fragment only large modules
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

DEGs_clusters <- split(top_gene_filtered$gene, top_gene_filtered$cluster)

# Function to count matches between DEGs and large modules
count_matches <- function(x, y) {
  sum(sapply(y, function(pattern) sum(x %in% pattern)))
}

# Apply the function to each element of DEGs_clusters using lapply
match_counts <- lapply(DEGs_clusters, function(sub_list) {
  sapply(large_modules, count_matches, y = sub_list)
})

# Convert the match_counts to a data frame
df <- as.data.frame(match_counts)

# Reshape the data frame to long format
df<-tibble::rownames_to_column(df,var = "Modules")
df_long <- tidyr::pivot_longer(df, cols = -Modules, 
                               names_to = "DEGs", values_to = "Coincidences")
# Calculating % based on modules
df_long <- df_long %>%
  mutate(Percentage = ifelse(Modules %in% names(large_modules), 
                             Coincidences / sapply(Modules, function(x) length(large_modules[[x]])) * 100,
                             NA)) %>%
  # Setting X and Y
  arrange(desc(Percentage))

# Setting axis order
fill_order <- sort(lengths(large_modules), decreasing = T)
df_long$Modules <- factor(df_long$Modules, levels = names(fill_order[c(1,3,2)]))

# Customizing labels
labels_vector <- setNames(paste0(names(fill_order)," (n=",as.numeric(fill_order),")"), 
                          names(fill_order))

ggplot(df_long, aes(x = DEGs, y = Percentage, fill = Modules)) +
  geom_bar(stat = "identity", color = "black") +  # Bars with black borders
  geom_text(aes(label = Coincidences), vjust = -0.5, size = 3) +  # Add Coincidences text above bars
  facet_wrap(~ Modules, scales = "free_x") +  # Facet wrap by Modules
  scale_fill_brewer(palette = "Set1",
                    # Use Set1 color palette for fill
                    labels = labels_vector) +  
  # Customize y-axis to show percentages with limits and breaks
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 40),  # Set y-axis limits
    breaks = seq(0, 40, by = 5)  # Set breaks at 5% intervals
  ) + 
  labs(
    title = "Coincidences between modules and DEGs per cluster",
    x = "DEGs per cluster",
    y = "Percentage of the module"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 14),  # Make facet titles larger
    panel.spacing = unit(1.5, "lines")  # Add spacing between facets
  )


#------------------------GENE ASSIGNATION
# Function to assign modules to clusters based on DEGs genes
gene_assignation <- function(DF_markers, selected_genes, 
                             log2fc_threshold = 0.6,
                             perform_second_round = FALSE,
                             n_rows = 2,
                             min_genes = 15) {
  library(dplyr)
  
  # Filter and group by gene
  top_gene_filtered <- DF_markers %>%
    filter(gene %in% selected_genes) %>%
    filter(avg_log2FC >= log2fc_threshold) %>%
    group_by(gene) %>%
    arrange(desc(avg_log2FC), desc(pct.1)) %>%
    slice(1:n_rows) %>%  # Select top n_rows per gene
    ungroup()
  
  # Optional second round of filtering
  if (perform_second_round) {
    top_gene_filtered <- top_gene_filtered %>%
      group_by(gene) %>%
      arrange(desc(avg_log2FC)) %>%
      mutate(diff_avg_log2FC = c(NA, diff(avg_log2FC))) %>%
      filter(abs(diff_avg_log2FC) <= 0.5 | is.na(diff_avg_log2FC)) %>%
      select(-diff_avg_log2FC) %>%  # Drop temporary column
      ungroup()
  }
  
  # Convert into a list of DEGs per cluster
  DEGs_list <- split(top_gene_filtered$gene, top_gene_filtered$cluster)
  
  # Keep only clusters with enough genes
  DEGs_list <- DEGs_list[lengths(DEGs_list) >= min_genes]
  
  return(DEGs_list)
}

#----------------------------WU DATASET: FRAGMENTING TURQUOISE MODULE
frag_turq_list<-gene_assignation(DF_markers = TAMDF,
                                 selected_genes = large_modules$turquoise,
                                 # min LFC to be considered Diff Express
                                 log2fc_threshold = 0.6,
                                 # Necessary if there is a difference of
                                 # more than 0.5 in fold change between clusters 
                                 perform_second_round = F,
                                 # Allow the gene to be classified in n clusters
                                 n_rows = 2,
                                 # Min group of genes to be considered a fragment
                                 min_genes = 15)
# Adding prefix
names(frag_turq_list)<-paste0("turquoise_",names(frag_turq_list))

# Visualize distribution of turquoise genes
library(venn)
frag_turq_inter <-attr(x = venn(frag_turq_list), "intersections")

# Turquoise is distributed between C2, C3, C5
# C3 and C5 share many genes while C5 has almost no unique DEGs.
# Merging C3 and C5 and compare against C2.
C3_C5vsC2 <-FindMarkers(TAMs,features = large_modules$turquoise,
                        ident.1 = c("C3","C5"), ident.2 = "C2",
                        # Calculate on normalize data
                        assay = "RNA", slot = "data", 
                        # Minimal fold change and expression %
                        logfc.threshold = 0.25, 
                        min.pct = 0.33)

# DEGs from both C3 and C5 constitute a fragment of larger module turquoise.
turq_C3_C5 <-C3_C5vsC2%>%
  filter(avg_log2FC>0&p_val_adj<0.01)%>%
  rownames()
# DEGs from C2 constitute a fragment of larger module turquoise.
turq_C2 <-C3_C5vsC2%>%
  filter(avg_log2FC<0&p_val_adj<0.01)%>%
  rownames()



#------------------------------FRAGMENTING BLUE MODULE
frag_blue_list<-gene_assignation(TAMDF, large_modules$blue,
                                 log2fc_threshold = 0.4,
                                 perform_second_round = F,
                                 n_rows = 2, min_genes = 15)
# Adding prefix
names(frag_blue_list)<-paste0("blue_",names(frag_blue_list))
# Visualize distribution of blue genes
library(venn)
frag_blue_inter<-attr(x = venn(frag_blue_list), "intersections")

# Blue is distributed between C0, C1, C2.
# C0 and C2 share many genes while C2 and C1 do not share genes.
# Merging C0 and C2 and compare against C1.
C0_C2vsC1<-FindMarkers(TAMs,features = large_modules$blue,
                       ident.1 = c("C0","C2"), ident.2 = "C1",
                       #Calculate on normalize data
                       assay = "RNA", slot = "data", 
                       #Minimal fold change and expression %
                       logfc.threshold = 0.25, #0.5?
                       min.pct = 0.33)

# DEGs from both C0 and C2 constitute a fragment of larger module blue.
blue_C0_C2 <-C0_C2vsC1%>%
  filter(avg_log2FC>0&p_val_adj<0.01)%>%
  rownames()
# DEGs from C1 constitute a fragment of larger module blue.
blue_C1 <-C0_C2vsC1%>%
  filter(avg_log2FC<0&p_val_adj<0.01)%>%
  rownames()


#------------------------------FRAGMENTING BROWN MODULE
# No gene from brown module is above 0.6 L2FC, unlike turquoise or blue
print(gene_assignation(TAMDF, large_modules$brown,
                       log2fc_threshold = 0.6,
                       perform_second_round = TRUE,
                       n_rows = 2, min_genes = 15))



#------------------------------DECISION
# Keep order from SortGenes by performing intersect function with Wu_SortGenes
fragments <-list(
  # Turquoise_C2
  intersect(Wu_SortGenes$turquoise,turq_C2),
  # turquoise_C3_C5
  intersect(Wu_SortGenes$turquoise,turq_C3_C5),
  # blue_C1
  intersect(Wu_SortGenes$blue,blue_C1),
  # blue_C0_C2
  intersect(Wu_SortGenes$blue,blue_C0_C2))

# Names
names(fragments) <-c("turquoise_C2","turquoise_C3&C5","blue_C1","blue_C0&C2")


# Discard brown module for being expressed in all clusters
Wu_fragmented <-Wu_SortGenes[!names(Wu_SortGenes)%in%names(large_modules)]

# Add fragmented large modules
Wu_fragmented<-c(Wu_fragmented,fragments)

# Change names
substitutions <- paste0("M",seq_along(Modules))
names(substitutions) <- paste0("^",names(sort(lengths(Modules), decreasing =T)))
names(Wu_fragmented) <- stringr::str_replace_all(names(Wu_fragmented), substitutions)

# saveRDS(Wu_fragmented, file = "D:/R_ScriptsPaper/Def_Objects/Wu_fragmented.rds")
#
# ============================================================================
# Section 2.2: FRAGMENTING LARGE MODULES WITH DEGs (AZIZI DATASET)
# ============================================================================
# AZIZI DATASET
Azizi_SortGenes <- readRDS("D:/R_ScriptsPaper/Def_Objects/Azizi_SortGenes.rds")
Modules <- Azizi_SortGenes

# Fragment only large modules
large_modules <- Modules[lengths(Modules) > 100]
# Marker of each cluster
TAMDF <-FindAllMarkers(TAMs,features = as.character(unlist(large_modules)),
                       # Calculate on normalize data
                       assay = "RNA", slot = "data", 
                       # Minimal fold change and expression %
                       logfc.threshold = 0.4, 
                       min.pct = 0.33, 
                       only.pos=T)

#----------------------------FRAGMENTING TURQUOISE MODULE
frag_turq_list <-gene_assignation(DF_markers = TAMDF,
                                  selected_genes = large_modules$turquoise,
                                  # min LFC to be considered Diff Express
                                  log2fc_threshold = 0.6,
                                  # Necessary if there is a difference of
                                  # more than 0.5 in fold change between clusters
                                  # It also means that there will be less genes 
                                  # assigned to 2 clusters
                                  perform_second_round = F,
                                  # Allow the gene to be classified in n clusters
                                  n_rows = 1,
                                  # Min group of genes to be considered a fragment
                                  min_genes = 15)

# Adding prefix
names(frag_turq_list) <-paste0("turquoise_",names(frag_turq_list))

# Visualize distribution of turquoise genes
library(venn)
frag_turq_inter <-attr(x = venn(frag_turq_list), "intersections")
#------------------------------DECISION
# Turquoise will be divided in two
turq_C5 <-frag_turq_list$turquoise_C5
#
turq_C4 <-frag_turq_list$turquoise_C4

# Keep order from SortGenes by performing intersect function with Wu_SortGenes
fragments <-list(
  # Turquoise_C5
  intersect(Azizi_SortGenes$turquoise,turq_C5),
  # turquoise_C4
  intersect(Azizi_SortGenes$turquoise,turq_C4)
  )

# Names
names(fragments) <-c("turquoise_C5","turquoise_C4")


# Discard brown module for being expressed in all clusters
Azizi_fragmented <-Azizi_SortGenes[!names(Azizi_SortGenes)%in%names(large_modules)]

# Add fragmented large modules
Azizi_fragmented <-c(Azizi_fragmented,fragments)

saveRDS(Azizi_fragmented, file = "D:/R_ScriptsPaper/Def_Objects/Azizi_fragmented.rds")
#
# ============================================================================
# Section 3: WHOLE MODULE ENRICHMENT AND DISTRIBUTION (BROWN)
# ============================================================================
# First module enrichment
library(GSVA)
datExpr<-GetAssayData(TAMs, assay = "RNA", slot = "data")
#Removing mitochondrial, Ribosomal and ORF genes
bad_patterns <- 'RP[0-9]|^MT|[0-9][0-9][0-9][0-9]|^RPL|^RPS|[0-9]orf[0-9]|^ATP'
a <- rownames(datExpr)[!grepl(bad_patterns, rownames(datExpr))]
datExpr<-datExpr[a,]
#Remove genes expressed in less than 5% of TAMs
datExpr<-datExpr[rowSums(datExpr == 0)/ncol(datExpr)<0.95,]
#Asses enrichment score with GSVA
TAM_enrich<-as.data.frame(gsva(datExpr,Modules,method="gsva",parallel.sz=1))

#Adding enrichment data to seurat object
library(dplyr)
ES_DF<-as.data.frame(t(TAM_enrich))
A<-colnames(ES_DF)
scores<-paste(A, "score", sep = "_")
colnames(ES_DF)<-scores

# CHOOSE RESPECTIVE DATASET

#saveRDS(ES_DF, file = "D:/R_ScriptsPaper/Def_Objects/Wu_ES_DF_WholeSig.rds")

#saveRDS(ES_DF, file = "D:/R_ScriptsPaper/Def_Objects/Azizi_ES_DF_WholeSig.rds")

#---- VISUALIZE DISTRIBUTION
#Import calculations for signatures
ES_DF <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_ES_DF_WholeSig.rds")

# Changing names
Modules <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_SortGenes.rds")
substitutions <- paste0("M",seq_along(Modules))
names(substitutions) <- paste0("^",names(sort(lengths(Modules), decreasing =T)))
colnames(ES_DF) <- stringr::str_replace_all(colnames(ES_DF), substitutions)

#Add to metadata
TAMs_scores <- AddMetaData(TAMs,metada=ES_DF)
# Only for
sc_selected <- as.character(substitutions)[grep("brown|blue|turquoise",names(substitutions))]
scores <- paste0(sc_selected,"_score")

# larger modules
scores <- c("M1_score","M2_score","M3_score")

# VIOLIN PLOT
# Function to calculate median point 
median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}
#Sorting by median in accordance with violin plots
obj.list <- SplitObject(TAMs_scores, split.by = "ident")

# Were to store violin plots
vln_list<-list()
#Loop for every score
for (ii in seq_along(scores)) {
  median_list<-vector()
  # One loop per cell type
  for (jj in seq_along(obj.list)){
    #Calculate median for each gene
    Q2<-as.numeric(quantile(obj.list[[jj]]@meta.data[,scores[ii]])["50%"])
    #Each member of the list contains a vector of every median for every gene 
    median_list[jj]<-Q2
  }
  #Convert list into dataframe
  names(median_list)<-names(obj.list)
  
  # Sort by median
  xOrder<-names(sort(median_list, decreasing = T))  
  #-----PLOT---
  vln_list[[ii]]<-VlnPlot(TAMs_scores, features = scores[ii], pt.size=0)+
    NoLegend() +
    # Order in X
    scale_x_discrete(limits=xOrder)+
    #labs(title = scores_aes[ii])+
    ylab("Enrichment score")+
    # More limit to add significant comparisons
    ylim(-1.1,1.3)+
    theme(axis.text.y = element_text(size = rel(0.7)),
          axis.text.x = element_text(size = rel(0.8)),
          #El texto horizontal
          axis.title.y = element_text(size = rel(0.8)),
          axis.title.x = element_blank(),
          plot.title = element_text(size = rel(0.8), family = "serif"),
          # Reduce space between plots
          plot.margin=unit(c(-0.5,-0.5,-0.5,0.5), "cm"))+
    # Add point at median
    stat_summary(fun = median.stat, geom='point', size = 1, colour = "black")
  print(paste("Loop",ii, sep = " "))
}

ggpubr::ggarrange(plotlist = vln_list, ncol=4, nrow = 2)



# BARR PLOT
# Calculate average ES in each cluster
obj.list <- SplitObject(TAMs_scores, split.by = "ident")
# Average ES for each score
mean.list<-list()
for (ii in seq_along(obj.list)) {
  mean.list[[ii]]<-apply(obj.list[[ii]]@meta.data[,scores], 2, mean)
}
names(mean.list)<-names(obj.list)

# Convert each element to a dataframe
mean_df <- do.call(rbind, mean.list)  
# Reshaping for plotting
long_mean <-reshape2::melt(mean_df)
colnames(long_mean)<- c("cluster","score","avg_ES")

# Ordering by cluster
long_mean$cluster<-factor(unique(long_mean$cluster), levels = paste0("C",c(0:7)))

ggplot(long_mean, aes(x = cluster, y = avg_ES, fill = score)) +
  geom_bar(stat = "identity", color = "black") +  # Bar plot with black borders
  facet_wrap(~ score) +        # Facet by Score
  labs(
    title = "Average score per cluster",
    x = "Seurat clusters",
    y = "Average ES"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1") +           # Use Set1 color palette for fill
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better visibility
    strip.text = element_text(size = 14),               # Larger facet titles
    legend.position = "none"                            # Remove legend (optional)
  )


# 
FeaturePlot(TAMs_scores, features = scores, ncol = 2,
            #max.cutoff = "q99",
            cols = rev(spect.pal))
#
# ============================================================================
# Section 4: DISCARD UNINTERESTING PATHWAYS (MKI67, CALR, ACTB)
# ============================================================================
Wu_fragmented <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_fragmented.rds")
# Add brown just in case
Wu_SortGenes <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_SortGenes.rds")
Wu_fragmented$M2 <- Wu_SortGenes$brown
names(Wu_fragmented) <-c("ISG15","FCN1","MKI67","CXCL9","CALR","FOLR2","APOE","IL1B","FOSB","ACTB")

# Load required libraries
library(org.Hs.eg.db)
library(limma)
library(dplyr)
library(tibble)
library(stringr)

# Download the list of pathways for KEGG
PathDF <- getKEGGPathwayNames(species = "hsa", remove.qualifier = TRUE)
colnames(PathDF) <- c("ID", "Pathway")

# Over representation analysis 
pathways_list <-list()
for (ii in seq_along(Wu_fragmented)) {
  Wu_as.df <- as.data.frame(Wu_fragmented[[ii]])
  Wu_as.df$EntrezID <- mapIds(org.Hs.eg.db, 
                           keys = Wu_as.df[, 1], 
                           column = "ENTREZID", 
                           keytype = "SYMBOL", 
                           multiVals = "first")
  
  # Enrichment with GO and KEGG databases
  goRest <- goana(Wu_as.df$EntrezID, species = "Hs")
  k <- kegga(Wu_as.df$EntrezID, species = "Hs")
  
  # Adding pathway names to KEGG
  k$ID <- str_remove_all(rownames(k), "path:")
  k$Pathway <- NULL
  k <- left_join(k, PathDF, by = "ID")
  k <- column_to_rownames(k, var = "ID")
  
  # Sorting pathways with p-value < 0.001 in GO and 0.01 in KEGG. More than 3 genes in path
  goSort <- goRest %>%
    filter(Ont == "BP" & P.DE < 10^-4 & DE >2) %>%
    arrange(P.DE)
  KEGGSort <- k %>%
    filter(P.DE < 0.01 & DE >2) %>%
    arrange(P.DE)
  KEGGSort$Ont <- rep("KEGG", nrow(KEGGSort))
  
  # Reordering columns of KEGGSort
  KEGGSort <- KEGGSort[, order(match(colnames(KEGGSort), c("Pathway", "Ont", "N", "DE", "P.DE")))]
  colnames(KEGGSort) <- colnames(goSort)
  # Binding results from both databases
  pathways_list[[ii]] <- rbind(goSort, KEGGSort)
  rownames_to_column(pathways_list[[ii]], var = "ID")
  print(paste0("Loop ",ii))
}

names(pathways_list)<-names(Wu_fragmented)

#---------------------GROUPING BY REDUNDANCY
library(stringdist)
library(dplyr)
library(stringr)

# Function to remove redundant terms with P.DE priority
remove_redundant_terms <- function(df, threshold, n.copies=2) {
  terms <- df$Term
  to_eliminate <- "positive |negative |regulation of |defense |response to | response$|signaling pathway| process$|cellular | cell$|cell |^trans"
  # Clean up the terms by removing redundant phrases
  terms<-str_remove_all(terms,to_eliminate)
  p_values <- df$P.DE
  
  # Compute a distance matrix based on string similarity
  dist_matrix <- stringdistmatrix(terms, terms, method = "jw")  # Jaro-Winkler similarity
  
  # Group terms based on similarity threshold
  # Higher threshold for more flexible clustering
  groups <- hclust(as.dist(dist_matrix)) %>%
    cutree(h = threshold)
  
  # Process dataframe to select representative terms
  result_df <- df %>%
    mutate(Group = groups) %>%
    add_count(Group, name = "PathwayCount") %>%
    group_by(Ont, Group) %>%
    arrange(P.DE, .by_group = TRUE) %>%  # Sort by P.DE within each group
    filter((n() >= 3 & row_number() <= n.copies) | (n() < 3 & row_number() == 1)) %>%
    ungroup()
  
  # Original names will be drawn
  return(result_df)
}
# Top ten per Data Base
filtered_terms <-lapply(pathways_list, remove_redundant_terms, 
                       threshold=0.3,
                       n.copies=1)

# Add the "Module" column with the corresponding name
for (name in names(filtered_terms)) {
  filtered_terms[[name]]$Module <- name
}

# Concatenating every DF into a single unified DF
go_kegg_DFresults <- do.call("rbind", filtered_terms)

# Save as csv
# write.csv(go_kegg_DFresults, file = "D:/R_ScriptsPaper/Def_Objects/Wu_go_kegg_DFresults.csv", row.names = F)

# Setting order
go_kegg_DFresults$Module <- factor(go_kegg_DFresults$Module ,
                                   levels = c("FCN1","ISG15","CXCL9","IL1B","FOSB","FOLR2","APOE","MKI67","CALR","ACTB"))

# Filter and choose number of pathways to represent 
go_kegg_filtered <- go_kegg_DFresults%>%
  # Group by Ont and globally arrange by P.DE and DE within Ont
  group_by(Ont, Module) %>%
  arrange(P.DE, desc(DE), .by_group = TRUE) %>%
  # Extract only the top pathways per subtype
  slice_head(n= 10) %>%
  ungroup()

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
          axis.ticks.x = element_blank()) +
    # Add vertical dashed lines
    geom_vline(xintercept = midpoints, linetype = "dashed")
  
}
# Plot side by side
go_kegg_plot <- cowplot::plot_grid(plotlist=PLIST, ncol = 2)
#
# VESTIGIOS####
#
# Further filtering / selecting pathways to represent for ... module
ISG15_paths<-c("reponse to virus","innate immune response",
               "response to type I interferon","response to cytokine",
               "positive regulation of interferon-beta production",
               "interleukin-27-mediated signaling pathway",
               "cytosolic pattern recognition receptor signaling pathway")

FCN1_paths<-c("innate immune response","response to stress","inflammatory response",
              "granulocyte migration","cell adhesion","response to cytokine")

MKI67_paths<-c("cell cycle process","chromosome organization","DNA metabolic process",
               "cellular component organization","DNA replication")

CXCL9_paths<-c("response to cytokine","response to type II interferon",
               "response to tumor necrosis factor","response to interleukin-1",
               "regulation of lymphocyte activation","neutrophil chemotaxis")

CALR_paths<-c("protein folding in endoplasmic reticulum","intracellular transport",
              "ERAD pathway","regulation of apoptotic process","Motor proteins")

APOE_paths<-c("antigen processing and presentation","complement activation",
              "regulation of hydrolase activity","response to lipoprotein particle",
              "cell junction disassembly","lipid metabolic process",
              # KEGG Paths
              "Lysosome","Complement and coagulation cascades",
              "Inflammatory bowel disease","Sphingolipid metabolism",
              "Cholesterol metabolism")

FOLR2_paths<-c("complement activation","cell junction disassembly","defense response",
               "Systemic lupus erythematosus"," Complement and coagulation cascades",
               "Efferocytosis","Phagosome")
C1<-c("response to lipopolysaccharide","apoptotic process","inflammatory response",
      "temperature homeostasis","neutrophil chemotaxis","endothelium development",
      "regulation of wound healing",
      "TNF signaling pathway","IL-17 signaling pathway","NF-kappa B signaling pathway")

FOSB_paths<-c("response to stress","response to unfolded protein",
              "response to temperature stimulus",
              "myeloid cell differentiation", "cellular response to interleukin-1",
              "eosinophil chemotaxis",
              "regulation of ERK1 and ERK2 cascade",
              "Estrogen signaling pathway","Th17 cell differentiation",
              "MAPK signaling pathway")
#
# ============================================================================
# Section 5: GENE FILTERING BASED ON MEMBERSHIP
# ============================================================================
# Apply only on Wu dataset because it is the only signatures that will be measured
Wu_fragmented <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_fragmented.rds")

# Function to select the top proportion of elements in each list item
top_proportion <- function(gene_list, proportion = 0.25) {
  lapply(gene_list, function(genes) {
    # Calculate the number of elements to retain based on the specified proportion
    n <- floor(proportion * length(genes))
    # Subset the first n elements (assumes they are already ordered)
    genes[1:n]
  })
}

# Reduce only modules with more than 30 genes
Wu_BestGenes <-top_proportion(Wu_fragmented[lengths(Wu_fragmented)>30], 
                              proportion = 0.4)
# Append smaller modules
Wu_BestGenes <-c(Wu_BestGenes,Wu_fragmented[lengths(Wu_fragmented)<=30])

names(Wu_BestGenes) <-c("ISG15","FCN1","MKI67","CXCL9","CALR","FOLR2","APOE","FOSB","IL1B")
# Setting desired order
Wu_BestGenes <-Wu_BestGenes[c("FCN1","ISG15","CXCL9","IL1B","FOSB","FOLR2","APOE")]

# saveRDS(Wu_BestGenes, file = "D:/R_ScriptsPaper/Def_Objects/Wu_BestGenes.rds")

#
