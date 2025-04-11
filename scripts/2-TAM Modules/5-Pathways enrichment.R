# ============================================================================
# Name: Fig1
# Author: elopez
# Date: Jun-30-2024
# Description: 
# TODO: 
# ============================================================================
library(Seurat)
library(dplyr)
# Analyzing only most polarized states
WuTAM_Pheno <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_Pheno.rds")
scores <-colnames(WuTAM_Pheno@meta.data)[grep("score",colnames(WuTAM_Pheno@meta.data))]
Sig.Name <-stringr::str_remove_all(scores, ".score")
# Selecting only the most polarized states
Polarized_TAMs <-subset(WuTAM_Pheno, idents=Sig.Name)

# Cluster order
levels(Polarized_TAMs) <- Sig.Name
# ============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================================
# DEA perform with normalized data (RNA, data)
# Both positive and negative FC to extract all relevant pathways.
# fgsea must be performed with a full list of foldchange (discarding non-significant paths)
DF_DEA<-FindAllMarkers(Polarized_TAMs, assay = "RNA",slot = "data", 
                       logfc.threshold = 0.1, min.pct = 0.1)

# Arrange by log2 fold change and only significant < 0.05
DEG_Poles<-DF_DEA %>%
  group_by(cluster)%>%
  arrange(desc(avg_log2FC), .by_group = T) %>% filter(p_val_adj < 0.05)

# Remove ribosomal an mithocondrial genes
bad_patterns <- 'RP[0-9]|^MT|[0-9][0-9][0-9][0-9]|^RPL|^RPS|[0-9]orf[0-9]|^ATP'
DEG_Poles<-DEG_Poles[!grepl(bad_patterns, DEG_Poles$gene),]
genes<-DEG_Poles$gene
#Convert gene name to ENTREZID
library("AnnotationDbi")
library("org.Hs.eg.db")
install.packages("C:/Users/erick/Downloads/org.Mm.eg.db_3.20.0.tar.gz", repos = NULL, type = "source")

DEG_Poles$EntrezID <- mapIds(org.Hs.eg.db,
                          keys=genes,
                          column="ENTREZID",
                          keytype="SYMBOL",
                          multiVals="first")

# Save DEA results for latter
# saveRDS(DEG_Poles, file = "D:/R_ScriptsPaper/Def_Objects/DEG_Poles.rds")
#
# ============================================================================
# GSEA
# ============================================================================
DEG_Poles <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/DEG_Poles.rds")
#Split genes by cluster
dea_list <-group_by(DEG_Poles, cluster, .add=T)%>% 
  group_split()
#Signature to use:
HGK <-readRDS(file = "D:/R_Scripts/EssentialObjects/HGK.rds")
library(fgsea)
ThreeDB<-list()
# Min and max pathway size for Hs.H and KEGG
for (i in 1:length(HGK)){
  #15 and 50 (stringent threshold) for GO pathways
  if (length(HGK[[i]])==7658){
    minS<-15
    maxS<-50
  }else{
    minS<-5
    maxS<-100
  }
  #Performing GSEA for every DEA gene from every module
  pathways <- HGK[[i]]
  topPathways<-list()
  for (G in 1:length(dea_list)){
    dfl<-dea_list[[G]]
    RANG <- dfl$avg_log2FC
    names(RANG) <- dfl$EntrezID
    #Remove NA in EntrezID
    RANG<-RANG[!is.na(names(RANG))]
    #Perform GSEA with pathways from Org.Hs
    prev_fgseaRes <- fgsea(pathways, RANG, minSize=minS, maxSize = maxS)
    #Only pathways significantly enriched
    fgseaRes<-filter(prev_fgseaRes, pval < 0.01)
    #top 10 overexpressed pathways
    topUp <- fgseaRes %>% 
      filter(ES > 0) %>% 
      top_n(10, wt=-pval)
    #top 10 underexpressed pathways
    topDown <- fgseaRes %>% 
      filter(ES < 0) %>% 
      top_n(10, wt=-pval)
    #Biding undder and overexpressed pathaways and sorting by NES
    topPathways[[G]] <- bind_rows(topUp, topDown) %>% 
      arrange(-NES)
  }
  names(topPathways) <-as.vector(unique(DEG_Poles$cluster))
  ThreeDB[[i]] <-topPathways
}
# Unifiying results from 3 databases
TOPPATHWAYS <-c(ThreeDB[[1]],ThreeDB[[2]],ThreeDB[[3]])

# Summary
TP <-ThreeDB[[2]]
for (i in 1:length(TP)){
  print(names(TP[i]))
  print(TP[[i]]$pathway)
  print(TP[[i]]$NES)
  print(TP[[i]]$size)
}
# saveRDS(TOPPATHWAYS, file = "D:/R_ScriptsPaper/Def_Objects/TOPPATHWAYS.rds")
#
# ============================================================================
# EXTRACTING LEADING EDGE FOR EACH PATHWAY
# ============================================================================
# Extracting leading edge of every pathways
TOPPATHWAYS <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/TOPPATHWAYS.rds")
Ledge<-list()
LE<-list()
for (P in 1:length(TOPPATHWAYS)){
  PN<-TOPPATHWAYS[[P]]$pathway
  a<-TOPPATHWAYS[[P]]$leadingEdge
  names(a)<-PN
  Ledge[[P]]<-a
  LE[[P]] <-PN
}
names(Ledge) <-names(TOPPATHWAYS)
# Simplifying leading edges duplicates
unif <-unlist(Ledge, recursive = FALSE, use.names = F)
names(unif) <-unlist(LE)
# Retrieving DE Genes
DEG_Poles <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/DEG_Poles.rds")
# Juntando vias en comun y convirtiendo ENTREZID en SYMBOL
Leading_edges <-split(unlist(unif, use.names = FALSE), rep(names(unif), lengths(unif)))
for (u in 1:length(Leading_edges)) {
  a<-unique(Leading_edges[[u]])
  aa<-as.data.frame(DEG_Poles[DEG_Poles$EntrezID%in%a,])
  Leading_edges[[u]]<-unique(aa$gene)
}
# saveRDS(Leading_edges,file = "D:/R_ScriptsPaper/Def_Objects/Leading_edges.rds")
#
# ============================================================================
# ENRICHMENT OF LEADING EDGE IN EACH POPULATION (PERFORMED ON CLUSTER)
# ============================================================================
# List of leading edge
Leading_edges <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Leading_edges.rds")
# Removing previous scores
scores <-colnames(WuTAM_Pheno@meta.data)[grep("score",colnames(WuTAM_Pheno@meta.data))]
WuTAM_Pheno@meta.data[,scores] <-NULL
datExpr <-GetAssayData(WuTAM_Pheno, assay = "RNA", slot = "data")
# Removing mitochondrial, Ribosomal and ORF genes
bad_patterns <- 'RP[0-9]|^MT|[0-9][0-9][0-9][0-9]|^RPL|^RPS|[0-9]orf[0-9]|^ATP'
a <- rownames(datExpr)[!grepl(bad_patterns, rownames(datExpr))]
datExpr<-datExpr[a,]
#Remove genes expressed in less than 5% of TAM_Mo
datExpr<-datExpr[rowSums(datExpr == 0)/ncol(datExpr)<0.95,]

library(GSVA)
TAM_enrich <-as.data.frame(gsva(datExpr,Leading_edges,method="gsva",parallel.sz=1))
ES_DF <-as.data.frame(t(TAM_enrich))

#saveRDS(ES_DF, file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_gsvaPaths.rds")
#
# ============================================================================
# REPRESENTATIVE PATHWAYS (GROUP BY REDUNDANCY)
# ============================================================================
GSVA_Res <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_gsvaPaths.rds")

#Adding enrichment data to seurat object
TAM_paths <-AddMetaData(WuTAM_Pheno, metada=GSVA_Res)

# In order to calculate average ES, we will split df by identity/Modules column
Idents(TAM_paths) <- TAM_paths$Modules

# Average ES in polarized pops
S.Obj <-subset(TAM_paths, idents =Sig.Name)

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

#Manually select representative pathways
# Define the patterns for each group
pattern_groups <- list(
  APOE = "ENDOCYTOSIS|LIPID_|LIPO|PPAR",
  FOLR2 = "COMPLEMENT|CHOND|CLATH|EPITHELIAL_CELL|GROWTH_FACTOR",
  CXCL9 = "INTERFERON|(?=.*ANTIGEN)(?=.*ENDOGENOUS)|STAT3|ALLOGRAFT",
  FCN1 = "ANGIOGENESIS|ANTIMICROBIAL",
  FOSB = "(?=.*CELLULAR)(?=.*STRESS)|HEAT|ERK|MAPK|MYC",
  ISG15 = "VIRAL",
  IL1B = "CELL_CHEMOTAXIS|KRAS|ERBB|TGF"
)

# Add a column indicating the group based on patterns
# Create a new column for group patterns
reduced_terms$GroupPattern <- "Unassigned"  # Default value

# Iterate through the pattern_groups and assign group names
for (group in names(pattern_groups)) {
  reduced_terms$GroupPattern[str_detect(reduced_terms$Term, pattern_groups[[group]])] <- group
}

# Hack to order GroupPattern acording to custom order
# Define the substitutions as a named character vector
# Names are the second values
substitutions <- c("FCN1" = "a_FCN1",
                   "ISG15" = "b_ISG15",
                   "CXCL9" = "c_CXCL9",
                   "IL1B" = "d_IL1B",
                   "FOSB" = "e_FOSB",
                   "FOLR2" = "f_FOLR2",
                   "APOE" = "g_APOE")

reduced_terms$GroupPattern <-stringr::str_replace_all(reduced_terms$GroupPattern, 
                                                       substitutions)

# Filter out "Unassigned" and all rows for ISG15 except the one with the highest value
CustomPaths <- reduced_terms %>%
  # Remove "Unassigned" rows
  filter(GroupPattern != "Unassigned") %>%
  group_by(GroupPattern) %>%
  # Keep only the top row for ISG15
  filter(!(GroupPattern == "b_ISG15" & HighestValue < max(HighestValue))) %>%
  arrange(desc(HighestValue)) %>%  # Ensure order by HighestValue within each group
  reframe(Term = Term)%>%
  pull(Term)

# saveRDS(CustomPaths, file = "D:/R_ScriptsPaper/Def_Objects/CustomPaths.rds")
