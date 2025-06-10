# ============================================================================
# Name: AcrossTissue Module assign
# Author: elopez
# Date: May-14-2024
# Description: 
# TODO: 
# ============================================================================
library(Seurat)
library(ggplot2)
library(dplyr)

#Seurat with estimated ES
Tissues_DF<-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Tissues_ES_DF.rds")
scores<-colnames(Tissues_DF)
Sig.Name<-stringr::str_remove_all(scores, ".score")

AcrossTissues<-readRDS(file="D:/R_ScriptsPaper/Objects/AcrossTissues_Wu.rds")
#
# ============================================================================
# Section 1: DETERMINING CELLS WITH HIGH ES
# ============================================================================
#Select cells with score>= x
#Select in two fractions for better characterization
#Select cells with ES higher than 0.x for each column
Highest<-apply(Tissues_DF,2,function(col) rownames(Tissues_DF)[col>=0.6])
names(Highest)<-Sig.Name

Snd.High<-apply(Tissues_DF,2,function(col) rownames(Tissues_DF)[col>=0.2 & col <0.6])
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
S.Obj_Pheno <-AddMetaData(AcrossTissues, Tissues_DF)
S.Obj_Pheno <-AddMetaData(S.Obj_Pheno, Phenotypes, col.name = "Modules")

#Cells with no high ES of any signature is Unassigned TAMs
S.Obj_Pheno$Modules<-ifelse(is.na(S.Obj_Pheno$Modules),
                            "Unassigned",
                            as.character(S.Obj_Pheno$Modules))
#
Idents(S.Obj_Pheno)<-S.Obj_Pheno@meta.data$Modules

# saveRDS(S.Obj_Pheno, file = "D:/R_ScriptsPaper/Def_Objects/Tissue_Pheno.rds")

#Saving elements for latter
#List of cells with high for each signature (No  phenotypes)
COMPLETE_ES_LIST<-list()
for (i in 1:length(Sig.Name)) {
  #Selecting phenotypes with each signature
  B<-Final.list[grep(Sig.Name[i],names(Final.list))]
  COMPLETE_ES_LIST[[i]]<-as.character(unlist(B))
}
names(COMPLETE_ES_LIST)<-Sig.Name

# saveRDS(COMPLETE_ES_LIST, file = "D:/R_ScriptsPaper/Def_Objects/TissueScores_HighCells_List.rds")
#