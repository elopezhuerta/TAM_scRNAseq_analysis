# =============================================================================
# Name: Module assigment
# Author: huerta
# Date: Apr-25-2024
# Description: Assignation of phenotypes based on quantile values.
# TODO:
# =============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)
#Load seurat object with calculated ES
#Wu DATASET
TAMs <-readRDS(file="D:/R_ScriptsPaper/Def_Objects/WuTAMs.rds")
GSVA_Res<-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_ES_DF_TAMs.rds")


#AZZI DATASET
TAMs <-readRDS(file="D:/R_ScriptsPaper/Def_Objects/Azizi_TAMs.rds")
GSVA_Res <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Azizi_ES_DF_TAM.rds")


#Add to metadata
Sig.Name <-stringr::str_remove_all(colnames(GSVA_Res), ".score")

#Palette
reds.pal<-RColorBrewer::brewer.pal(9,"YlOrRd")
spect.pal<-rev(RColorBrewer::brewer.pal(11,"Spectral"))
greens.pal<-RColorBrewer::brewer.pal(8,'Greens')
#
# ============================================================================
# Section 1: DETERMINING CELLS WITH HIGH ES
# ============================================================================
#Select cells with score>= x
#Select in two fractions for better characterization
#Select cells with ES higher than 0.x for each column
Highest<-apply(GSVA_Res,2,function(col) rownames(GSVA_Res)[col>=0.6])
names(Highest)<-Sig.Name

Snd.High<-apply(GSVA_Res,2,function(col) rownames(GSVA_Res)[col>=0.2 & col <0.6])
names(Snd.High)<-Sig.Name
# ============================================================================
# Section 2: Defining intersection in both list
# =============================================================================
#Finding cells with highest ES for more than one signature with intersections
library(UpSetR)
#Observe modules with no intersection
complete.interDF <-upset(fromList(Highest),nsets = length(Highest))
#Obtanin intersections as df of binary data
binar.modules <-complete.interDF$New_data
#Setting rownames
rownames(binar.modules) <-unique((unlist(Highest, use.names=F)))
#Selecting columns with 1
mod.inter <-vector()
for (i in 1:nrow(binar.modules)) {
  s<-binar.modules[i,]==1
  mod.inter[i]<-paste(colnames(binar.modules)[s], collapse = ":")
}
#Module combined with clusters
HIGH_inter<-split(rownames(binar.modules),mod.inter)

#PERFOM SAME ANALYSIS WIT CELLS WITH THE SECOND HIGHEST VALUE
#Some cells with the highest values might be repeated in cells with second highest values for a different signature
#Remove repeated cells in Highest and Second highest group (Snd.High)
Snd.unique  <-lapply(Snd.High, setdiff, unlist(Highest))

#Finding cells with second highest ES for more than one signature with intersections
#Observe modules with no intersection
complete.interDF<-upset(fromList(Snd.unique),nsets = length(Snd.unique))
#Obtanin intersections as df of binary data
binar.modules<-complete.interDF$New_data
#Setting rownames
rownames(binar.modules)<-unique((unlist(Snd.unique, use.names=F)))
#Selecting columns with 1
mod.inter<-vector()
for (i in 1:nrow(binar.modules)) {
  s<-binar.modules[i,]==1
  mod.inter[i]<-paste(colnames(binar.modules)[s], collapse = ":")
}
#Module combined with clusters
SECOND_inter<-split(rownames(binar.modules),mod.inter)


#JOINING BOTH LISTS
#Determine phenotypes unique to each list
uniq.pheno<-c(setdiff(names(HIGH_inter), names(SECOND_inter)),
              setdiff(names(SECOND_inter), names(HIGH_inter)))

#Save phenotypes unique to each list and removing NULL elements
uniq.list<-purrr::discard(c(HIGH_inter[uniq.pheno],SECOND_inter[uniq.pheno]), 
                          is.null)
## IN ORDER TO UNIFY HIGH_inter and SECOND_inter, BOTH LIST MOST HAVE THE SAME ELEMTS
## THEREFORE, WE MUST REMOVE UNIQUE PTHENOTYPES IN EITHER LISTS
## UNIQUE PHENOTYPES WILL BE APPENDE AFTER UNIFICATION OF CLEANED HIGH_ AND SECOND_ inter
HIGH_inter<-HIGH_inter[!names(HIGH_inter)%in%uniq.pheno]
SECOND_inter<-SECOND_inter[!names(SECOND_inter)%in%uniq.pheno]
#Joining list
unified_list<-Map(c, HIGH_inter, SECOND_inter)
# Recovering list of unique phenotypes
Final.list<-c(unified_list, uniq.list)

#To add this element to seurat object, it needs to be a named vector
named.list<-list()
for (i in c(1:length(Final.list))){
  A<-names(Final.list[i])
  Cell.IDs<-Final.list[[i]]
  #Name each cell with the respective name
  B<-rep(A,length(Cell.IDs))
  names(B)<-Cell.IDs
  named.list[[i]]<-B
}
Phenotypes<-unlist(named.list)

# Add to metadata
S.Obj_Pheno <-AddMetaData(TAMs, metada=GSVA_Res)
S.Obj_Pheno <-AddMetaData(S.Obj_Pheno, Phenotypes, col.name = "Modules")

#Cells with no high ES of any signature is Unassigned TAMs
S.Obj_Pheno$Modules<-ifelse(is.na(S.Obj_Pheno$Modules),
                            "Unassigned",
                            as.character(S.Obj_Pheno$Modules))
#
Idents(S.Obj_Pheno)<-S.Obj_Pheno@meta.data$Modules

# CHOOSE LOCATION
# WU DATASET
# saveRDS(S.Obj_Pheno, file = "D:/R_ScriptsPaper/Def_Objects/WuTAM_Pheno.rds")

# Azizi DATASET
# saveRDS(S.Obj_Pheno, file = "D:/R_ScriptsPaper/Def_Objects/AziziTAM_Pheno.rds")

#Saving elements for latter
#List of cells with high for each signature (No  phenotypes)
COMPLETE_ES_LIST<-list()
for (i in 1:length(Sig.Name)) {
  #Selecting phenotypes with each signature
  B<-Final.list[grep(Sig.Name[i],names(Final.list))]
  COMPLETE_ES_LIST[[i]]<-as.character(unlist(B))
}
names(COMPLETE_ES_LIST)<-Sig.Name

#WU DATASET
# saveRDS(COMPLETE_ES_LIST, file = "D:/R_ScriptsPaper/Def_Objects/WuScores_HighCells_List.rds")

#AZIZI DATASET
# saveRDS(COMPLETE_ES_LIST, file = "D:/R_ScriptsPaper/Def_Objects/AziziScores_HighCells_List.rds")
# ============================================================================
# Section 3: CORRELATION PLOTS (FAILED)
# =============================================================================
plot_correlation <- function(high_cells, ES_DF, cell_type1, cell_type2) {
  
  # Extract cell names for the two selected cell types
  cells1 <- high_cells[[cell_type1]]
  cells2 <- high_cells[[cell_type2]]
  
  Col_1<-paste0(cell_type1, "_score")
  Col_2<-paste0(cell_type2, "_score")
  
  # Subset ES_DF for cells in both cell types
  correlation_df <- ES_DF[rownames(ES_DF) %in% cells1 | rownames(ES_DF) %in% cells2, c(Col_1,Col_2)]
  
  colnames(correlation_df)<-c("CellType1_score","CellType2_score")
  
  # Create and return the scatter plot
  ggplot(correlation_df, aes(x = CellType1_score, y = CellType2_score)) +
    geom_point(alpha = 0.7, color = "blue") +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    labs(
      title = paste("Correlation between", cell_type1, "and", cell_type2, "Scores"),
      x = paste(cell_type1, "Score"),
      y = paste(cell_type2, "Score")
    ) +
    theme_minimal()
}
library(ggplot2)
plot_correlation(ES_DF= GSVA_Res,high_cells = COMPLETE_ES_LIST,
                 cell_type1 = "JUN", cell_type2 = "FOSB")
#