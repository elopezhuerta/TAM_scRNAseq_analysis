# ============================================================================
# Name: Fig2-Enrichment in TME
# Author: elopez
# Date: Aug-04-2024
# Description: 
# TODO: 
# ============================================================================
library(Seurat)
# ES in TME
TME_scores <-readRDS(file = "D:/Cluster_Scripts/Descargar/GSVA/TestWuTME_DF_scores.rds")
# Extracting scores and polarized states names
scores <-colnames(TME_scores)
single.modules <-stringr::str_remove_all(scores, ".score")
TAM_turns <-paste0("TAM_", single.modules)

# Aesthetic
scores_aes <- stringr::str_replace_all(scores,"_"," ")

# Import and edit object with TILs
Infilt_scores <-readRDS(file = "D:/R_ScriptsPaper/Objects/WuBRCA_Refined.rds")
Infilt_scores <-AddMetaData(Infilt_scores, metadata = TME_scores)



# Loading object with TAM identity
TAM_Modules <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_Pheno.rds")
# Only most polarized states
Poles <- subset(TAM_Modules,idents = single.modules)

# Replacing TAM identities while pasting TAM prefix
Infilt_scores<-SetIdent(Infilt_scores, cells = names(Poles$Modules), 
                   value = paste0("TAM_",as.character(Poles$Modules)))

# Removing remaining TAMs and Doublets
Infilt_scores<-subset(Infilt_scores, idents = c("TAMs","Cycling TAMs", "Doublets"), 
                      invert=T)
# Identities should be the same as celltypes
Infilt_scores$celltypes<-Idents(Infilt_scores)
#
# ============================================================================
# SHAPIRO TEST (NORMALITY TEST)
# ============================================================================
# Test by score per cell type
TME_by_celltypes <-dplyr::group_split(Infilt_scores@meta.data[,c("celltypes",scores)], celltypes)

# Create a df indicating normality for each cell type per score
for (ii in seq_along(TME_by_celltypes)) {
  print(ii)
  DFL<-TME_by_celltypes[[ii]]
  #If cell type is > 5k, make sample because 5k is limit for shapiro test
  if (nrow(DFL)>5000) {
    DFL <- DFL[sample(nrow(DFL), 5000), ]
  }
  #Test for each module signature, while removing celltype column
  shapiro_results <- apply(DFL[,-1], 2, shapiro.test)
  # Pval for each score
  shap_res<-vector()
  for (jj in 1:ncol(DFL[,-1])) {
    # If pval<0.05, data is not normally distributed
    shap_res[jj]<-shapiro_results[[jj]]$p.value>0.05
  }
  result_df<-data.frame(names(shapiro_results),shap_res)
  colnames(result_df)<-c(as.character(unique(DFL$celltypes)),"Normal_distribution")
  print(result_df)
}
#
# ============================================================================
# KRUSKAL-WALLIS
# ============================================================================
library(dplyr)
library(ggpubr)
library(FSA)  # For Dunn's Test
#Populations of interest
# Were to store Dunn results
ns.comparisons <- list()

for (ii in seq_along(scores)) {
    # Only TAM population of interest
    S.Obj <- subset(Infilt_scores,idents =TAM_turns[-ii], invert=TRUE)
    # Create a df with only the desired column
    df <- S.Obj@meta.data[,c("celltypes",scores[ii])]
    # Homogenize names
    colnames(df)<-c("celltypes","ES")
    # Perform the Kruskal-Wallis Test
    kruskal_test <- kruskal.test(ES ~ celltypes, data = df)
    
    # If the Kruskal-Wallis test is significant, perform Dunn's Test for post-hoc analysis
    if (kruskal_test$p.value < 0.05) {
      dunn_test <- dunnTest(ES ~ celltypes, data = df, method = "bonferroni")
      dunn_res <- dunn_test$res
      # Rows were single.modules[ii] appear only once
      #logic_vector<-grepl(single.modules[ii], dunn_res$Comparison)
      logic_vector <- grepl("TAM", dunn_res$Comparison)
      # Display only ns differences involving the respective TAM population
      ns.pair <- dunn_res$Comparison[logic_vector&dunn_res$P.adj > 0.01]
      if (length(ns.pair)==0) {
        ns.pair<-"All significant"
      }
      ns.comparisons[[ii]] <- ns.pair
    } else {
      print(paste0("Kruskal-Wallis test was n.s.; no post-hoc test is performed for ",scores[ii]))
    }
  }
# Summary of n.s. comparison with the respective population (single.modules[ii])
names(ns.comparisons) <- scores
#
# ============================================================================
# VIOLIN PLOTS
# ============================================================================
library(ggplot2)
#Funcion para marcar punto medio del violin
median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}
# Were to store violin plots
vln_list <- list()

# Loop for every score (TAM_turns coincide) to sort by median values
for (ii in seq_along(scores)) {
  median_list <- vector()
  # Perform the operations with all TAMs except one
  S.Obj <- subset(Infilt_scores,idents =TAM_turns[-ii], invert=TRUE)
  #Sorting by median in accordance with violin plots
  obj.list <- SplitObject(S.Obj, split.by = "ident")
  # One loop per cell type
  for (jj in seq_along(obj.list)){
    #Calculate median for each gene
    Q2<-as.numeric(quantile(obj.list[[jj]]@meta.data[,scores[ii]])["50%"])
    #Each member of the list contains a vector of every median for every gene 
    median_list[jj]  <- Q2
  }
  #Convert list into dataframe
  names(median_list) <- names(obj.list)
  # Sort by median
  xOrder <- names(sort(median_list, decreasing = T))
  
  #-----PLOT---
  vln_list[[ii]] <- VlnPlot(Infilt_scores, features = scores[ii], pt.size=0)+
    NoLegend() +
    # Order in X
    scale_x_discrete(limits=xOrder)+
    labs(title = scores_aes[ii])+
    ylab("ES")+
    # More limit to add significant comparisons
    ylim(-1.1,1.3)+
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 6),
          #El texto horizontal
          axis.title.y = element_text(size = 8),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 10, face = "plain"),
          # Reduce space between plots
          plot.margin=unit(c(-0.5,-0.5,-0.5,0.5), "cm"))+
    # Add point at median
    stat_summary(fun = median.stat, geom='point', size = 1, colour = "black")
  print(paste("Loop",ii, sep = " "))
}
# Plot array
jpeg("D:/R_ScriptsPaper/Figures/Manuscrito/pre_Fig5.png", 
     res = 600,width = 18, height = 17.25, units = "cm")
ggpubr::ggarrange(plotlist = vln_list, ncol=3, nrow = 3)
dev.off()

# Display significant comparisons
print(ns.comparisons)
#
# ============================================================================
# TAM SEPARATION
# ============================================================================
library(dplyr)
library(ggpubr)
library(FSA)  # For Dunn's Test
#Populations of interest
# Were to store Dunn results
ns.TAM_comparisons <- list()
# Only TAM population of interest
TAM.Obj <- subset(Infilt_scores,idents =TAM_turns)

for (ii in seq_along(scores)) {
  # Create a df with only the desired column (celltypes are idents)
  df <-TAM.Obj@meta.data[,c("celltypes",scores[ii])]
  # Homogenize names
  colnames(df) <-c("celltypes","ES")
  # Perform the Kruskal-Wallis Test
  kruskal_test <- kruskal.test(ES ~ celltypes, data = df)
  
  # If the Kruskal-Wallis test is significant, perform Dunn's Test for post-hoc analysis
  if (kruskal_test$p.value < 0.05) {
    dunn_test <- dunnTest(ES ~ celltypes, data = df, method = "bonferroni")
    dunn_res <-dunn_test$res
    # Rows were single.modules[ii] appear only once
    logic_vector <-grepl(single.modules[ii], dunn_res$Comparison)
    #logic_vector<-grepl("TAM", dunn_res$Comparison)
    # Display only ns differences involving the respective TAM population
    ns.pair <-dunn_res$Comparison[logic_vector&dunn_res$P.adj>0.01]
    if (length(ns.pair)==0) {
      ns.pair <- "All significant"
    }
    ns.TAM_comparisons[[ii]]<-ns.pair
  } else {
    print(paste0("Kruskal-Wallis test was n.s.; no post-hoc test is performed for ",scores[ii]))
  }
}

# Summary of n.s. comparison with the respective population (single.modules[ii])
names(ns.TAM_comparisons) <- scores
#

library(ggplot2)
#Funcion para marcar punto medio del violin
median.stat <- function(x){
  out <- quantile(x, probs = c(0.5))
  names(out) <- c("ymed")
  return(out) 
}
#Sorting by median in accordance with violin plots
obj.list <- SplitObject(TAM.Obj, split.by = "ident")
# Were to store violin plots
TAM_vln_list <- list()

#Loop for every score
for (ii in seq_along(scores)) {
  median_list<-vector()
  #S.Obj <- subset(Infilt_scores,idents =TAM_turns[-ii], invert=TRUE)
  #Sorting by median in accordance with violin plots
  #obj.list <- SplitObject(S.Obj, split.by = "ident")
  # One loop per cell type
  for (jj in seq_along(obj.list)){
    #Calculate median for each gene
    Q2 <- as.numeric(quantile(obj.list[[jj]]@meta.data[,scores[ii]])["50%"])
    #Each member of the list contains a vector of every median for every gene 
    median_list[jj] <- Q2
  }
  #Convert list into dataframe
  names(median_list) <- names(obj.list)
  # Sort by median
  xOrder <- names(sort(median_list, decreasing = T))
  
  #-----PLOT---
  TAM_vln_list[[ii]]<-VlnPlot(Infilt_scores, features = scores[ii], pt.size=0)+
    NoLegend() +
    # Order in X
    scale_x_discrete(limits=xOrder)+
    labs(title = scores_aes[ii])+
    ylab("ES")+
    # More limit to add significant comparisons
    ylim(-1.1,1.3)+
    theme(axis.text.y = element_text(size = rel(0.7)),
          axis.text.x = element_text(size = rel(0.6)),
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
# Plot array
ggpubr::ggarrange(plotlist = TAM_vln_list, ncol=4, nrow = 2)
# Display significant comparisons
print(ns.TAM_comparisons)
#
# ============================================================================
# Joining
# ============================================================================
upper_part <-ggpubr::ggarrange(plotlist = vln_list, ncol=4, nrow = 2)

lower_section <-ggpubr::ggarrange(plotlist = TAM_vln_list, 
                                  ncol=length(TAM_vln_list), nrow = 1)

cowplot::plot_grid(upper_part,lower_section, nrow = 2,
                   labels = c("a","b"), label_size = 20,
                   rel_heights = c(1,0.5))
