#setwd("E:/Maestr?a/cluster_juriquilla/out/QC")
library(Seurat)
#Batch correction#####
QC_merg<-readRDS("QC_merg.rds")
Merg_list <- SplitObject(QC_merg, split.by = "orig.ident")
#CHECAR METODOS DE NORMALIZACION?
for (i in 1:length(Merg_list)) {
  Merg_list[[i]] <- NormalizeData(Merg_list[[i]], verbose = FALSE)
  Merg_list[[i]] <- FindVariableFeatures(Merg_list[[i]], selection.method = "vst", 
                                         nfeatures = 2000, verbose = FALSE)
}
#QUE ME DIJO DE BATCH CORRECTION? CHECAR METODOS TAMBIEN
#k.filter es 50 porque el menor (BC08_N) es de 59 c?lulas
Merg.anchors <- FindIntegrationAnchors(object.list = Merg_list, dims = 1:30,k.filter=50)
#IntegrateData returns new assay with 'batch-corrected' expression matrix for all cells, enabling them to be jointly analyzed.
Merg.anchors <- IntegrateData(anchorset = Merg.anchors, dims = 1:30)
# switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(Merg.anchors) <- "integrated"
#QUE MAS SE HACE CON PCA??
Merg.anchors <- ScaleData(Merg.anchors, verbose = FALSE)
#DEBUGGING
na_cells <- colSums(as.matrix(is.na(Merg.anchors@assays$integrated@data)))
to_remove <- names(na_cells)[na_cells>0]
if(length(to_remove) > 0){
  print("Found some cells with NA values. They will be removed. These are the cells:")
  print(length(to_remove))
  Merg.anchors@meta.data$keep <- TRUE
  Merg.anchors@meta.data[to_remove, "keep"] <- FALSE
  Merg.anchors <- subset(Merg.anchors, keep)
  Merg.anchors@meta.data <- Merg.anchors@meta.data[, -c(which(colnames(Merg.anchors@meta.data) == "keep"))]
}
#Necessary to re-run ScaleData
Merg.anchors <- ScaleData(Merg.anchors, verbose = FALSE)
Merg.anchors <- RunPCA(Merg.anchors, npcs = 30, verbose = FALSE)
#NO A?ADIR TSNE PORQUE SE NECESITA UNO NUEVO PARA TAMS
Merg.anchors <- RunUMAP(Merg.anchors, reduction = "pca", dims = 1:30)
Merg.anchors <- FindNeighbors(Merg.anchors, dims = 1:10)
Merg.anchors <- FindClusters(Merg.anchors, resolution = 0.5)

#Salvar por si hay que repetir anotacion
saveRDS(Merg.anchors, file="Merg_anchors.rds")

#Annotation####
#De la librer?a celldex obtuvimos la referencia con la funcion:
#hpca.se <- HumanPrimaryCellAtlasData()
hpca.se<-readRDS(file = "hpca.rds")
#ES LA NORMALIZADA???
hESCs<-Merg.anchors[["RNA"]]@data
#Quitando genes poco interesantes o que estan definidos como ORF 
bad_patterns <- 'RP[0-9]|^MT|[0-9][0-9][0-9][0-9]|^RPL|^RPS'
a <- rownames(hESCs)[!grepl(bad_patterns, rownames(hESCs))]
hESCs<-hESCs[a,]

library(SingleR)
#Anotando por celula
pred.hesc <- SingleR(test = hESCs, ref = hpca.se,
                     method = "single",
                     assay.type.test=1,
                     labels = hpca.se$label.main)

#Celulas cuya anotaci?n no coincidi? con al menos otras 49 celulas (eventos aislados) ser?n unknown
library(dplyr)
unknown<-as.data.frame(table(pred.hesc$pruned.labels))
unknown<-as.character(filter(unknown,Freq<=100)$Var1)
#Cambiando anotaci?n
new.labels<-ifelse(pred.hesc$pruned.labels%in%unknown,"unknown",pred.hesc$pruned.labels)
pred.hesc$labels<-new.labels

#Para saber todas las labels antes de reclasificar como unknown
saveRDS(pred.hesc, file = "pred_hesc.rds")

#Para saber que c?lula corresponde a cada tipo celular.
Clust.annot<-as.data.frame(cbind(rownames(pred.hesc),pred.hesc$labels))
#Assigning cell type lables to cells
new.cluster.ids <- Clust.annot$V2
names(new.cluster.ids) <- Clust.annot$V1
Idents(Merg.anchors)<-new.cluster.ids
#Guardando celltypes en metadata
Merg.anchors$celltypes<-Idents(Merg.anchors)
#Drop low quality annotated cells
all.cells<-names(table(Idents(Merg.anchors)))
Annotated<-subset(Merg.anchors,idents = all.cells[!all.cells%in%"unknown"])

saveRDS(Annotated, file="Annotated.rds")
