# =============================================================================
# Name: Preprocessing BRCA Azzi files
# Author: huerta
# Date: Apr-17-2024
# Description: Read tsv files and convert them into seurat objects
# TODO: Run untar function only once to avoid duplicated files.
# =============================================================================

#Trabajar en la ruta /mnt/Genoma/drobles/elopez/scRNAseq
ruta_final<-"D:R_ScriptsPaper/Raw data"
#Correr la siguiente funcion "untar()" una sola vez!!, si no se duplican los archivos
#untar("D:/Cluster_Juriquilla/in/GSE114725_RAW.tar",exdir =ruta_final)
#Para tener en una variable todos los archivos que se extrajeron de GSE114725_RAW.tar
archivos_gz<-list.files(path = ruta_final)
#Combinar el nombre del archivo con su respectiva ruta
muestras<-paste(ruta_final,archivos_gz,sep = "/")
#Seleccionar muestras tumorales y control
muest_select<-c("NORMAL", "TUMOR", "BLOOD")
#Todo lo que es normal o tumor (16 muestras)
descargar_estas<-grep(paste(muest_select,collapse = "|"),muestras,value = T)
#install.packages("data.table")
#install.packages("R.utils")
library(R.utils)
library(data.table)
#La siguiente funcion (pre_treat) abarca lo siguiente: 
#cambiar filas por columnas (celulas en columnas y genes en filas)
#asignar celulas como colnames y eliminar fila de cells duplicada
pre_treat<-function(i){
  mydataframe<-fread(descargar_estas[[i]],header = T)
  mydataframe[is.na(mydataframe)] <- 0
  tdf<-t(mydataframe)
  colnames(tdf) <- tdf[1, ]
  tdf<-tdf[-1,]
}
#La funci?n anterior se repetir? las veces que le indiquemos 1:length(descargar_estas) 
#y generar? una lista de matrices modificadas (DF_list) como le indicamos con la funci?n pre_treat
DF_list<-lapply(1:length(descargar_estas), FUN = "pre_treat")
#Preparando el nombre de project para cada objeto Seurat
library(stringr)
#quitar extensi?n de archivo
proj_name<-str_remove(descargar_estas, "_counts.csv.gz")
#Quitando ruta + el identificador /GSMXXXXXXX_
proj_name<-substr(proj_name,nchar(ruta_final)+13,nchar(proj_name))

library(Seurat)
#La funci?n Z crea el objeto Seurat como se especifica 
#(deja fuera celulas con menos de 200 genes y que eso no se repita en a menos 2 celulas mas)
z<-function(i){
  my_seurat<-CreateSeuratObject(DF_list[[i]], project = proj_name[i],
                                min.cells = 3, min.features = 200)
}
#Obtenemos una lista de n objetos Seurat (uno por muestra)
my_seurat<-lapply(1:length(DF_list),FUN="z")
saveRDS(my_seurat, file = "D:/R_Scripts/EssentialObjects/AzziRaw.rds")
