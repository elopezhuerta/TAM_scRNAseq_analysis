library(Seurat)
#Este objeto seurat fue creado con script Archivos comprimidos. Consta de las muestras BC06, BC07, BC08 BC03(BC03 consta de un tumor y un ctrl)
my_seurat<-readRDS("my_seurat.rds")
#Setting identifier based on sample origin
proj_name<-vector()
for (i in 1:length(my_seurat)) {
  proj_name[i]<-unique(as.character(my_seurat[[i]]$orig.ident))
}
halb_seurat<-my_seurat[c(-1)]

#Se a?ade peque?o identificador al ID de cada celula para que no este repetido (mismo batch)
merg_seurat<-merge(x=my_seurat[[1]], y=halb_seurat, add.cell.ids = c(proj_name[1], proj_name[c(-1)]))
#La muestra de donde proviene esta en orig.ident
#QC
#Se calcula %mitocondriales
merg_seurat[["percent.mt"]] <- PercentageFeatureSet(merg_seurat, pattern = "^MT-")
#Identificando muestras con cantidades inusualmente altas o bajas (3nmads):
#Como genes detectados (nFeature), library size (nCount), % MT(da?adas)
library(scater)
#library(scuttle)
lots_genes<-isOutlier(merg_seurat@meta.data$nFeature_RNA,type = "higher",
                      batch =merg_seurat$orig.ident)
damaged<-isOutlier(merg_seurat@meta.data$percent.mt,type="higher",
                   batch =merg_seurat$orig.ident)
lots_counts<-isOutlier(merg_seurat@meta.data$nCount_RNA,type="higher",
                          batch =merg_seurat$orig.ident)
#Quitando outliers
QC_merg<-merg_seurat[,!(lots_genes|damaged|lots_counts)]

# Add in a metadata column that indicates sample/condition
library(stringr)
#Para pre-QC
neuen<-as.character(Idents(merg_seurat))
sampID<-substr(neuen,1,6)
merg_seurat$Samp_Cond <- sampID
Idents(merg_seurat)<-sampID
saveRDS(merg_seurat, file = "merg_seurat.rds")
#Post-QC
neuen_QC<-as.character(Idents(QC_merg))
sampID_QC<-substr(neuen_QC,1,6)
QC_merg$Samp_Cond <- sampID_QC
Idents(QC_merg)<-sampID_QC
saveRDS(QC_merg, file = "QC_merg.rds")

#Diagnostico visual
library(ggplot2)
library(cowplot)
library(patchwork)
#Que no haya dobletes o c?lulas muy desviadas
pre<-FeatureScatter(merg_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+
  labs(title = "pre-QC")
post<-FeatureScatter(QC_merg, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ 
  labs(title = "post-QC")
#Set single common legend
both <- pre + post & theme(legend.position = "bottom")
jpeg("duplets.jpg", res = 300,width = 465,height = 225, units = "mm")
both + plot_layout(guides = "collect")
dev.off()

#Hay que corroborar que el pico de los violinplots sea menos largo
#Funcion para hacer violines acomodados en filas (stacked violin seurat)####
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(0.9), angle = 90), 
          axis.text.y = element_text(size = rel(0.8)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(size = rel(0.9), angle = 45), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

#Continua####
#Salvando
jpeg("pre_QC.jpg", res = 300,width = 465,height = 225, units = "mm")
StackedVlnPlot(merg_seurat, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
dev.off()
jpeg("post_QC.jpg", res = 300,width = 465,height = 225, units = "mm")
StackedVlnPlot(QC_merg, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
dev.off()

#Informacion basica
a<-as.matrix(table(QC_merg$orig.ident)/table(merg_seurat$orig.ident)*100)
b<-as.matrix(table(QC_merg$orig.ident))
c<-as.matrix(table(Idents(QC_merg)))
ab<-cbind(a,b[,1])
n<-max(length(ab),length(c))
length(ab) <- n                      
length(c) <- n
abc<-as.data.frame(cbind(ab,c),optional = T)
colnames(abc)<-c("%kept_cells","num total cell","total samp/cond")
write.csv(abc, file = "basic_info.csv")
