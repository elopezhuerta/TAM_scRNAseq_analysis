# ============================================================================
# SUPPLEMENTARY TABLE S1: MODULE GENES
# ============================================================================
#-----------SHEET 1: ORIGINAL MODULES
#Signature from WuWGCNA sorted by KME
WuModules <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_SortGenes.rds")
Wu_datKME <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_datKME.rds")

#Signature from Azizi sorted by KME
AziziModules <- readRDS("D:/R_ScriptsPaper/Def_Objects/Azizi_SortGenes.rds")
Azizi_datKME <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Azizi_datKME.rds")


Module_Excel <- function(Modules, datKME, dataset){
  # Initialize an empty list to store results
  module_dfs <- list()
  # Loop through each module color
  for (color in names(Modules)) {
    # Get genes for current module
    module_genes <- Modules[[color]]
    color_modify <- stringr::str_remove(color,"anti_") 
    # Subset datKME for these genes and the current module's column
    subset_df <- datKME[module_genes, color_modify, drop = FALSE]
    # Add module color column
    subset_df$Original_module_name <- color
    # Rename the Membership_value column consistently
    colnames(subset_df)[1] <- "Membership_value"
    subset_df <- tibble::rownames_to_column(subset_df,"Gene")
    subset_df$Dataset <- dataset
    # Store in list
    module_dfs[[color]] <- subset_df
  }
  # Combine all modules into one dataframe
  My_Excel <- dplyr::bind_rows(module_dfs)
}
Wu_Excel<-Module_Excel(WuModules, Wu_datKME, dataset = "Wu, 2021")
Azizi_Excel<-Module_Excel(AziziModules, Azizi_datKME, dataset = "Azizi, 2018")

# 
TableS1_Sheet1 <-rbind(Wu_Excel,Azizi_Excel)

#-----------SHEET 2: REFINED MODULES
Wu_fragmented <- readRDS("D:/R_ScriptsPaper/Def_Objects/Wu_fragmented.rds")
names(Wu_fragmented) <- c("ISG15","FCN1","MKI67","CXCL9","CALR","FOLR2",
                          "APOE","IL1B","FOSB")
# As DF
Wu_fragDF <- plyr::ldply(Wu_fragmented,cbind)
colnames(Wu_fragDF) <- c("New_name","Gene")

Wu_refined<- dplyr::left_join(Wu_fragDF, Wu_Excel, by="Gene")


Azizi_fragmented <- readRDS("D:/R_ScriptsPaper/Def_Objects/Azizi_fragmented.rds")
names(Azizi_fragmented) <- c("FOLR2","APOE","HLA_II","FCN1","ISG15")
# As DF
Azizi_fragDF <- plyr::ldply(Azizi_fragmented,cbind)
colnames(Azizi_fragDF) <- c("New_name","Gene")
# Adding the rest of the data
Azizi_refined<- dplyr::left_join(Azizi_fragDF, Azizi_Excel, by="Gene")

# Into a single sheet
TableS1_Sheet2 <-rbind(Wu_refined,Azizi_refined)

#-----------SHEET 3: MODULE SIGNATURES
Wu_BestGenes <- readRDS("D:/R_ScriptsPaper/Def_Objects/Wu_BestGenes.rds")

TableS1_Sheet3 <- plyr::ldply(Wu_BestGenes,cbind)
colnames(TableS1_Sheet3)<-c("Module signature name","Signature")
#-----------SHEET 4: CLUSTER-ASSOCIATED DEGs
library(Seurat)
#Wu DATASET
Wu_TAMs <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAMs.rds")
# Setting clusters as idents
Idents(Wu_TAMs) <-paste0("C",Wu_TAMs$seurat_clusters)
# Cluster associated DEGs
TableS1_Sheet4 <-FindAllMarkers(Wu_TAMs,
                       # Calculate on normalize data
                       assay = "RNA", slot = "data", 
                       # Minimal fold change and expression %
                       logfc.threshold = 0.4, 
                       min.pct = 0.33, 
                       only.pos=T)

TableS1_Sheet4 <- TableS1_Sheet4[TableS1_Sheet4$p_val_adj<0.01,]

# AZIZI DATASET
Azizi_TAMs <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Azizi_TAMs.rds")
# Setting clusters as idents
Idents(Azizi_TAMs) <-paste0("C",Azizi_TAMs$seurat_clusters)
# Cluster associated DEGs
TableS1_Sheet5 <-FindAllMarkers(Azizi_TAMs,
                                # Calculate on normalize data
                                assay = "RNA", slot = "data", 
                                # Minimal fold change and expression %
                                logfc.threshold = 0.4, 
                                min.pct = 0.33, 
                                only.pos=T)

TableS1_Sheet5 <- TableS1_Sheet5[TableS1_Sheet5$p_val_adj<0.01,]

#-----------UNIFYING ALL SHEETS
library(writexl)
# Create a named list of dataframes
excel_sheets <- list(
  "WGCNA derived modules" = TableS1_Sheet1,
  "Refined modules" = TableS1_Sheet2,
  "Module signatures" = TableS1_Sheet3,
  "Wu_Cluster associated DEGs" = TableS1_Sheet4,
  "Azizi_Cluster associated DEGs" = TableS1_Sheet5
)
# Write to Excel
write_xlsx(excel_sheets, "D:/R_ScriptsPaper/Tablas/SupplementaryTableS1.xlsx")
#
# ============================================================================
# SUPPLEMENTARY TABLE S2: MODULE PATHWAYS AND SIGNATURES
# ============================================================================
###-----------SHEET 1: OVERREPRESENTED PATHWAYS
Over_repr <- read.csv(file = "D:/R_ScriptsPaper/Def_Objects/Wu_go_kegg_DFresults.csv")
Over_repr_sheet <- Over_repr[,c(1:5,8)]

colnames(Over_repr_sheet) <- c("Overrepresented pathway","Database","Total genes","Gene hits",
                     "pval adjusted","Wu refined modules")
TableS2_Sheet1 <- Over_repr_sheet


###-----------SHEET 2: DIFF_EXPR POLARIZED TAMs 
DEG_Poles <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/DEG_Poles.rds")
TableS2_Sheet2 <- as.data.frame(DEG_Poles[,c(6,7,2,5,3,4)])
colnames(TableS2_Sheet2) <- c("TAM subset","Gene","avg_log2FC","pval_adj",
                             "%_express_respective_TAM","%_express_other_TAMs")
# Only positive markers
TableS2_Sheet2 <- TableS2_Sheet2[TableS2_Sheet2$avg_log2FC>0,]


##-----------SHEET 3: TAM SUBSETS ENRICHED PATHWAYS
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
#CHOOSE DATASET
#Wu, 2021 Dataset
WuTAM_Pheno <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_Pheno.rds")

scores <-colnames(WuTAM_Pheno@meta.data)[grep("score",colnames(WuTAM_Pheno@meta.data))]
Sig.Name <-stringr::str_remove_all(scores, ".score")

# 157 pathways GSVA
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
promsPaths <- plyr::ldply(mean_ES,rbind)
# We need TAM subset column as rownames to transprose table
TableS2_Sheet3 <- tibble::column_to_rownames(promsPaths, var = ".id")
TableS2_Sheet3 <- as.data.frame(t(TableS2_Sheet3))
TableS2_Sheet3 <- tibble::rownames_to_column(TableS2_Sheet3)
# 
colnames(TableS2_Sheet3) <- c("Enriched pathway", paste0(colnames(TableS2_Sheet3)[-1]," (average enrichment score)"))

###-----------SHEET 4: PATHWAY-ASSOCIATED SIGNATURES
Leading_edges <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Leading_edges.rds")
# Transform list in dataframe
TableS2_Sheet4 <- plyr::ldply(Leading_edges, cbind)
colnames(TableS2_Sheet4) <- c("Pathway", "Leading edge genes")

# Create a named list of dataframes
TableS2_excel_sheets <- list(
  "Overrepresented pathways" = TableS2_Sheet1,
  "Diff_expr of Polarized TAMs" = TableS2_Sheet2,
  "TAM subsets enriched pathways" = TableS2_Sheet3,
  "Pathway associated signatures" = TableS2_Sheet4
)
# Write to Excel
writexl::write_xlsx(TableS2_excel_sheets, "D:/R_ScriptsPaper/Tablas/SupplementaryTableS2.xlsx")
#
# ============================================================================
# SUPPLEMENTARY TABLE S3: PATHWAY-ASSOCIATED SIGNATURES
# ============================================================================
###-----------SHEET 1: INFLAMMATORY/TISSUE REPAIR SIGNATURES
UNI <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Leading_edges.rds")

# Tissue repair
repair_patterns <- "WOUNDING|ENDOCYTOSIS|DIFFERENTIATION|COAGULATION|LYSOSOME"
Repair <- unique(unlist(UNI[grep(repair_patterns,names(UNI))]))
# Interferon high
IFN_Alfa <- unique(unlist(UNI[grep("VIRAL|ALPHA|TYPE_I",names(UNI))]))
IFN_Gamma <- unique(unlist(UNI[grep("GAMMA|TYPE_II",names(UNI))]))
# Immune activation
Immune <- unique(unlist(UNI[grep("HUMORAL|KAPPAB|INFLAM",names(UNI))]))

# Unifying list
Resum <- list(Immune, IFN_Alfa, IFN_Gamma, Repair)
names(Resum) <- c("Innate_immunity","IFN_I","IFN_II","Tissue_repair")

# Creating df and Excel sheet
TableS3_Sheet1 <- plyr::ldply(Resum, cbind)
colnames(TableS3_Sheet1) <- c("Infl. or Repair signature", "Gene")
#
###-----------SHEET 2: Curated M1 M2 markers
# Excel is in Table directory. Add Manually

###-----------SHEET 3: Fold change values M1 M2
#Importing signature to test
file_path<-"D:/R_ScriptsPaper/Tablas/M1M2.xlsx"
df_M1M2 <- readxl::read_excel(file_path,col_names = T, range = "A1:B58")
colnames(df_M1M2) <- c("Marker","gene")

DEG_Poles <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/DEG_Poles.rds")
DEG_Poles <- as.data.frame(DEG_Poles[,c(6,7,2,5,3,4)])

# Filter and arrange the data frame
DEG_Poles <- DEG_Poles %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  filter(avg_log2FC > 0, # Only positive avg_log2FC values
         # Filter by genes in M1_M2
         #gene %in% M1_M2_genes,
         # Filter by adjusted p-value < 0.05
         p_val_adj < 0.05)


TableS3_Sheet3 <- dplyr::left_join(df_M1M2, DEG_Poles, by="gene")

colnames(TableS3_Sheet3) <- c("M1_M2 Marker","Gene","TAM subset","avg_log2FC","pval_adj",
                              "%_express_respective_TAM","%_express_other_TAMs")

# Replace NA in non-numeric columns with "n.s."
TableS3_Sheet3[] <- lapply(TableS3_Sheet3, function(x) {
  if (!is.numeric(x)) {
    x[is.na(x)] <- "n.s."
  }
  return(x)
})


###-----------SHEET 4: Tissue comparisons pvals
DF_list <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/complete_tissue_comparisons.rds")
TableS3_Sheet4 <- do.call("rbind",DF_list)
TableS3_Sheet4 <- tibble::rownames_to_column(TableS3_Sheet4, var = "Signature_score")
# Remove numeric suffixes
TableS3_Sheet4$Signature_score <- gsub("\\.\\d+$", "", TableS3_Sheet4$Signature_score)
# Change _ for vs
TableS3_Sheet4$Comparison <- gsub(" - "," vs ", TableS3_Sheet4$Comparison)

# Create a named list of dataframes
TableS3_excel_sheets <- list(
  "Infl_Tissue Repair Signatures" = TableS3_Sheet1,
  "FC values M1 M2" = TableS3_Sheet3,
  "Tissue comparisons" = TableS3_Sheet4
)
# Write to Excel
writexl::write_xlsx(TableS3_excel_sheets, "D:/R_ScriptsPaper/Tablas/SupplementaryTableS3.xlsx")
#
# ============================================================================
# SUPPLEMENTARY TABLE S4: FOSB DEGs
# ============================================================================

# DEA results are in  D:/R_ScriptsPaper/Def_Objects/DEG_FOSBvsAll.csv
#