# =============================================================================
# Name: Supplementary Figure S7 
# Author: elopez
# Date: 28-08-2024
# Description: 
# TODO: 
# ============================================================================
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
# Abundance estimated by Bayes Prism
BP.theta <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/BP_Wu_Poles.rds")

TAM_pops <- c("FCN1","ISG15","CXCL9","IL1B","FOSB","FOLR2","APOE")

# % format
BP.theta <- BP.theta*100

# Determining populations of interest
TAM_stringpattern <- grep(paste(TAM_pops, collapse = "|"), 
                          colnames(BP.theta), value = TRUE)

#Subsetting TAM results (% scale)
TAM_estimation <- BP.theta[,TAM_stringpattern]

#Subsetting TAM with log10 scale of 100%
TAM_log_scale <- log10(BP.theta[,TAM_stringpattern])

# CLEANING clinical dataframe
# Import clinical data
clinical_ready <- readRDS(file = "D:/R_Scripts/EssentialObjects/clinical_ready.rds")
# Choosing clinical parameter to analyze
All_Param <- c("PAM50_Subtype","Pathologic_stage","PR_Status","ER_Status","Her2_Status")

# Check variables to clean
apply(clinical_ready[,All_Param], 2, table)

# Cleaning columns from clinical data 
to_clean <- c("NA","Normal","[Discrepancy]","[Not Available]","[Not Evaluated]","Indeterminate","Equivocal")

for (i in seq_along(to_clean)) {
  clinical_ready[All_Param][clinical_ready[All_Param]==to_clean[i]] <- NA
}
#
# ============================================================================
# Fig. S7A: BOXPLOT ALL POPULATION ABUNDANCES IN THE TME 
# ============================================================================
# Calculate the xth percentile for each column in BP.theta
top_90th <- apply(BP.theta, 2, function(x) quantile(x, 0.9, na.rm = TRUE))

# Sort the variables based on the 75th percentile in decreasing order
sorted_top_90th <- names(top_90th)[order(top_90th, decreasing = TRUE)]

# Reshape the dataframe to long format
df_long <- reshape2::melt(BP.theta)

# Plot with the x-axis ordered by the 75th percentile values
TME_estim_boxes <- ggplot(df_long, aes(x = Var2, y = value, fill = Var2)) +
  geom_boxplot() +
  theme_minimal() +
  # Order x axis by 75th percentile
  scale_x_discrete(limits = sorted_top_90th) +
  labs(title = "Estimated cell type abundance in the tumor") +
  ylab("Percentage") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 13, face = "plain"),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  ) +  # Remove grid lines
  scale_fill_brewer(palette = "Set3")
#
# ============================================================================
# Fig. S7B-D: TAM ABUNDANCES WITHIN PATIENTS
# ============================================================================
# ABUNDANCE (%) OF TAM POPS
Only_TAM_abun <- as.data.frame(TAM_estimation)
Only_TAM_abun$Total_TAMs <- apply(Only_TAM_abun,1,sum)

# Calculate the 75th percentile (3rd quartile) for each column in BP.theta
TAM_75th <- apply(Only_TAM_abun, 2, function(x) quantile(x, 0.9, na.rm = TRUE))
# Sort the variables based on the 75th percentile in decreasing order
sorted_TAM_75th <- names(TAM_75th)[order(TAM_75th, decreasing = TRUE)]

# BOXPLOT ABUNDANT TAMs
gg_Only_TAM_abun <- reshape2::melt(Only_TAM_abun)

TAM_boxplots <- ggplot(gg_Only_TAM_abun, aes(x = variable, y = value, fill = variable)) +  # Add 'fill' to color the boxes
  geom_boxplot() +
  theme_minimal() +
  scale_x_discrete(limits = sorted_TAM_75th) +
  labs(title="Estimated TAM subset abundance in the tumor")+
  ylab("Percentage") + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 13, face = "plain"),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  ) +  # Remove grid lines
  scale_fill_brewer(palette = "Set3")


# CO-EXISTANCE WITHIN PATIENTS
library(pheatmap)
library(gplots)
#Palette color
spect.pal <- RColorBrewer::brewer.pal(11, "Spectral")

# Only top 4 abundant TAM populations
TAM_est_DF <- TAM_estimation[,sorted_TAM_75th[2:5]]
# Heatmap of %
TAM_coex_patients <- pheatmap(t(TAM_est_DF), col=rev(spect.pal), 
                              cluster_cols=T, cluster_rows=T, 
                            fontsize_row = 10,
                            border_color=NA, 
                            show_colnames = F,
                            annotation_legend = F, 
                            annotation_names_row = F)

# As ggplot object to arrange with other plots
TAM_coex_patients <- ggplotify::as.ggplot(TAM_coex_patients)

# CO-EXISTANCE PAIRS BULK
# CHOOSE APPROPIATE LIMIT FOR EACH POPULATION OF INTEREST
apply(TAM_est_DF, 2, quantile,c(.75,.80,.90,.95,1))
#Select cells with score>= x
#Select in two fractions for better characterization
#Select cells with ES higher than 0.x for each column
Highest<-apply(TAM_est_DF,2,function(col) rownames(TAM_est_DF)[col>=3])
names(Highest)<-TAMs_surv

Snd.High<-apply(TAM_est_DF,2,function(col) rownames(TAM_est_DF)[col>=1 & col <3])
names(Snd.High)<-TAMs_surv

# Adjusting for APOE
Snd.High$TAM_APOE<-rownames(TAM_est_DF)[TAM_est_DF[,"TAM_APOE"]>=0.5 & TAM_est_DF[,"TAM_APOE"] <3]

# Create the AllCombinations matrix
All.Combinations <- combn(names(Highest), 2)

#Frequency of module co-existance
createDotFreqPlot <- function(Modules, AllCombinations, dataset) {
  Frequency <- vector()
  Pairs <- vector()
  for (ii in 1:ncol(AllCombinations)) {
    # Frequency of high patients for the respective pair
    Highest_counts <- intersect(Highest[[AllCombinations[1,ii]]],
                                Highest[[AllCombinations[2,ii]]])
    
    Snd_counts <- intersect(Snd.High[[AllCombinations[1,ii]]],
                            Snd.High[[AllCombinations[2,ii]]])
    Frequency[ii] <- length(Highest_counts)+length(Snd_counts)
    # Combination name
    pre_pairs <- str_remove_all(All.Combinations[1:2, ii],"TAM_")
    Pairs[ii] <- paste(pre_pairs, collapse = " & ")
  }
  #Result as DF for plotting
  CombDF <- data.frame(Pairs, Frequency)
  CombDF$Pairs<-stringr::str_replace_all(CombDF$Pairs,"_","-")
  # Order by frequency
  CombDF <- CombDF[order(CombDF$Frequency), ]
  # Fix order for plotting by setting factor and levels
  CombDF$Pairs <- factor(CombDF$Pairs, levels = unique(CombDF$Pairs))
  
  DotFreq <- ggplot(CombDF, aes(x = Frequency, y = Pairs)) +
    geom_point(colour = "red", size = 2) +
    expand_limits(x = max(Frequency) + 15) +
    ylab("Co-existing TAMs in patients") +
    xlab("Number of High-patients") +
    geom_text(label = CombDF$Frequency, hjust = -1, size = 3.5) +
    theme_classic() +
    theme(axis.title = element_text(size = 13),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          # Adjust size and center
          plot.title = element_blank(),
          panel.border = element_rect(fill = NA)) +
    ggtitle(title)
  
  return(DotFreq)
}

patient_pair <-createDotFreqPlot(Highest, AllCombinations = All.Combinations, dataset = "(TCGA)")
#
# ============================================================================
# SURVIVAL: DETERMINING CELLS WITH HIGH ESTIMATED ABUNDANCE
# ============================================================================
# Only abundant TAMs
TAMs_surv <- c("TAM_CXCL9","TAM_ISG15","TAM_FOLR2","TAM_APOE")
#-------REMOVE NA ROWS BEFORE PERFORMING INTERSECTION
patient_surv <- clinical_ready$patient[!is.na(clinical_ready$Evento)]

df_to_clean <- as.data.frame(TAM_estimation[rownames(TAM_estimation)%in%patient_surv,TAMs_surv])

# CHOOSE APPROPIATE LIMIT FOR EACH POPULATION OF INTEREST
apply(df_to_clean, 2, quantile,c(.75,.80,.90,.95,1))
#Select cells with score>= x
#Select in two fractions for better characterization
#Select cells with ES higher than 0.x for each column
Highest<-apply(df_to_clean,2,function(col) rownames(df_to_clean)[col>=4])
names(Highest)<-TAMs_surv

Snd.High<-apply(df_to_clean,2,function(col) rownames(df_to_clean)[col>=2 & col <4])
names(Snd.High)<-TAMs_surv


Highest$TAM_APOE<-rownames(df_to_clean)[df_to_clean$TAM_APOE>=4]

Snd.High$TAM_APOE<-rownames(df_to_clean)[df_to_clean$TAM_APOE>=0.5 & df_to_clean$TAM_APOE <4]

#--------Defining intersection in both list
#Finding cells with highest ES for more than one signature with intersections
library(UpSetR)
#Observe modules with no intersection
complete.interDF<-upset(fromList(Highest),nsets = length(Highest))
#Obtanin intersections as df of binary data
binar.modules<-complete.interDF$New_data
#Setting rownames
rownames(binar.modules)<-unique((unlist(Highest, use.names=F)))
#Selecting columns with 1
mod.inter<-vector()
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
#Remove unique phenotypes from original lists
HIGH_inter<-HIGH_inter[!names(HIGH_inter)%in%uniq.pheno]
SECOND_inter<-SECOND_inter[!names(SECOND_inter)%in%uniq.pheno]
#Joining list
unified_list<-Map(c, HIGH_inter, SECOND_inter)
#Adding unique phenotypes to fina list
Final.list<-c(unified_list, uniq.list)
#
# ============================================================================
# Fig. S7E: SURVIVAL CURVES
# ============================================================================
# Selecting patients with high abundance of the phenotype of interest
my_list <- Final.list[TAMs_surv[!TAMs_surv%in%"TAM_CXCL9"]]

#----------SURVIVAL
# Convert the list to a named vector
abundance_subtitution <- unlist(my_list)
# Adjust the names
names(abundance_subtitution) <- c(rep(names(my_list), lengths(my_list)))

#name vector as dataframe
abund_df<-data.frame(patient= as.vector(abundance_subtitution),
                     Abundance= names(abundance_subtitution))
# Joining df
clinical_merge_abundance<-merge(clinical_ready,abund_df,
                                by='patient',all=TRUE)

# Only necessary columns
curves_DF<-clinical_merge_abundance[,c("Months","Evento","Abundance")]

# Removing NA from 
curves_DF<-curves_DF[!is.na(curves_DF$Evento)&!is.na(curves_DF$Abundance),]

curves_DF$Abundance <- factor(curves_DF$Abundance, levels = c("TAM_ISG15", 
                                                              "TAM_APOE", 
                                                              "TAM_FOLR2"))

#Survival analisis
a.km <- survfit(Surv(as.numeric(Months), Evento) ~ Abundance, data = curves_DF, type = "kaplan-meier")
#Para obtener median survival
#surv_median(a.km)
TAM_Curves <- ggsurvplot(fit = a.km, data = curves_DF, 
                       conf.int = FALSE, pval = F, 
                       surv.median.line = "none",
                       break.x.by = 50,
                       risk.table = "nrisk_cumevents",
                       risk.table.y.text.col = TRUE, 
                       risk.table.y.text = FALSE,
                       # Adjust this to make the numbers smaller
                       risk.table.fontsize = 4,  
                       legend.title = element_blank(),
                       legend.labs = c("High in ISG15 TAMs",
                                       "High in APOE TAMs", 
                                       "High in FOLR2 TAMs"),
                       # Position the legend to the right
                       legend = "right",  
                       # Use theme_classic() to remove the grid
                       ggtheme = theme_classic() +
                         theme(panel.grid = element_blank(),
                               legend.text = element_text(size = 10),
                               axis.title = element_text(size = 13),
                               axis.text = element_text(size = 10),
                               # Keep x and y axis lines
                               axis.line = element_line(color = "black"))
)

# Extract the main plot and the risk table
Surv_plot <- TAM_Curves$plot
Surv_table <- TAM_Curves$table +
  theme(plot.title = element_text(size = 12))

# Combine the survival plot and the risk table into one plot
combined_surv_plot <- cowplot::plot_grid(Surv_plot, Surv_table, 
                                         ncol = 1, rel_heights = c(2, 1))

# Perform pairwise log-rank tests to get pairwise p-values
pairwise_pvals <- pairwise_survdiff(Surv(as.numeric(Months), Evento) ~ Abundance, data = curves_DF)

# Extract the p-value matrix
print(pairwise_pvals$p.value)
#
# ============================================================================
# Fig. S7F: VOLCANO PLOT FOSB vs TOTAL MONOCYTES
# ============================================================================
# Import DEA results
DEGs_DF <- read.csv(file = "D:/R_ScriptsPaper/Def_Objects/DEG_FOSBvsAll.csv",
                    header = TRUE, row.names = 1)

#Intersecting signatures + rescued modules
BestGenes<-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Wu_BestGenes.rds")


genes_to_display<- c("FOSB","FOS","DUSP1","CCL3","NFKBIA","IL1B","SRGN","CXCL8",
                     "CHMP1B","CXCR4","ZFP36","RGS2")

#c("ID2","CD52","CEBPB","JUND")
#sign_genes <- head(intersect(rownames(DEGs_DF),BestGenes$FOSB), 10)
#no_sign_genes <- head(setdiff(rownames(DEGs_DF),BestGenes$FOSB),10)
# must be a dataframe
#genes_to_display <-unique(c(sign_genes, no_sign_genes))

# must be a dataframe
top_genes<-DEGs_DF[rownames(DEGs_DF)%in%genes_to_display,]

# Create a basic volcano plot
volcano_plot <- ggplot(DEGs_DF, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = ifelse(p_val_adj < 0.01, "Significant", "Not Significant")), 
             size = 0.5) +
  #Puntos rojos si es significativo
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
  theme_minimal() +
  labs(x = "Log2FC", y = "-log10(p-value)") +
  #Linea que indica significancia
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0)+
  theme(legend.position = "none", # Remove the legend
        panel.grid=element_blank(),
        # Add margins (top, right, bottom, left)
        plot.margin = unit(c(0.5, 0.5, 0.2, 0.2), "cm"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.line.x = element_line(color = "black",linewidth = 0.5),  # Customize x-axis line color and thickness
  )  
# No overlap in genes
library(ggrepel)
vcplot <- volcano_plot +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)),
                  size = 4,  # Adjust the size value as needed
                  box.padding = 0.3, point.padding = 0.1)

# Adding blank space to violin plot
space_vcplot <-cowplot::plot_grid(vcplot, NULL, nrow=2 ,rel_heights = c(1,0.2))
#
# ============================================================================
# JOINING
# ============================================================================
# Middle panel
midle <- cowplot::plot_grid(TAM_boxplots,TAM_coex_patients,patient_pair,
                          ncol = 3,
                          rel_widths = c(0.6,1,0.8))

Lower_Panel <- cowplot::plot_grid(combined_surv_plot, space_vcplot,ncol = 2)

# Final arrangement
cowplot::plot_grid(TME_estim_boxes, midle, Lower_Panel, nrow = 3,
                   rel_heights = c(0.5,0.6,0.8))

#
# ============================================================================
# LEFT OVERS: CLINICAL ASSOCIATIONS-BOXPLOT
# ============================================================================
# Merge the data frames, using row names of TAM_log_scale as the patient column
clin_TAMlog <- merge(clinical_ready[,c("patient",All_Param,"Evento")], TAM_log_scale, 
                     by.x = "patient", by.y = "row.names")

# Eliminating rows with no estimated abundance for TAM populations
long_clinTAM <- clin_TAMlog %>%
  pivot_longer(cols = TAM_HSPs:TAM_CXCL9,    # Specify the range of columns to pivot
               names_to = "variable",
               values_to = "value") 

# Setting TAM population order
long_clinTAM$variable<-factor(long_clinTAM$variable,
                              levels = sorted_TAM_75th[-1])
# Setting order of PAM50 subtype
long_clinTAM$PAM50_Subtype<-factor(long_clinTAM$PAM50_Subtype,
                                   levels = c("LumA","LumB","Her2","Basal"))

boxpl_list<-list()
for (ii in seq_along(All_Param)) {
  df.ggp<-long_clinTAM[,c(All_Param[ii],"variable","value")]
  colnames(df.ggp)<-c("Parameter","TAM_pop","Abundance")
  # Removing NA rows in each case
  df.ggp<-df.ggp[!is.na(df.ggp$Parameter),]
  boxpl_list[[ii]]<-ggboxplot(df.ggp, x = "TAM_pop", y = "Abundance",
                              fill = "Parameter",color = "Black")+
    #Reset background to white
    theme_classic()+
    ylab("% of tumor composition (log10)")+
    labs(title = All_Param[ii])+
    ylim(min(df.ggp$Abundance),max(df.ggp$Abundance)+0.8)+
    theme(legend.title = element_blank(),
          # Panel style
          panel.border = element_rect(fill=NA),
          axis.title = element_text(face = "bold",size = 11,family = "serif"),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1,vjust = 1),
          plot.title = element_text(face = "bold",size = 13,hjust = 0.5, 
                                    family = "serif", colour = "red"))
}

Clin_Asso_Plot <- cowplot::plot_grid(boxpl_list[[1]],boxpl_list[[2]],
                                     labels = c("d","e"), label_size = 15, 
                                     ncol = 2,rel_widths = 1)

# ----------STATISTICAL TESTS: KRUSKAL-WALLIS
# Merge the data frames, using row names of TAM_estimation as the patient column
clin_TAMestim <- merge(clinical_ready[,c("patient",All_Param,"Evento")], TAM_estimation, 
                       by.x = "patient", by.y = "row.names")

library(FSA)  # For Dunn's Test
Stage_PAM<-c("PAM50_Subtype","Pathologic_stage")
# Were to store Dunn results
DF_list<-list()
for (ii in seq_along(Stage_PAM)) {
  ns.comparisons<-list()
  for (jj in seq_along(TAM_pops)) {
    # Create a df with only the desired column
    sub_df<-clin_TAMestim[,c(Stage_PAM[ii],TAM_pops[jj])]
    colnames(sub_df)<-c("Parameter","Abundance")
    df_kruskal<-sub_df[!is.na(sub_df$Parameter), ]
    # Perform the Kruskal-Wallis Test
    kruskal_test <- kruskal.test(Abundance ~ Parameter, data = df_kruskal)
    # If the Kruskal-Wallis test is significant, perform Dunn's Test for post-hoc analysis
    if (kruskal_test$p.value < 0.01) {
      dunn_test <- dunnTest(Abundance ~ Parameter, data = df_kruskal, method = "bonferroni")
      dunn_res<-dunn_test$res
      # Display only ns differences involving the respective TAM population
      ns.pair<-dunn_res$Comparison[dunn_res$P.adj<0.01]
      ns.comparisons[[jj]]<-ns.pair
    } else {
      print(paste0("Kruskal-Wallis test was n.s.; no post-hoc test is performed for ",
                   Stage_PAM[ii], ": ",TAM_pops[jj]))
      ns.comparisons[[jj]]<-"n.s."
    }
  }
  names(ns.comparisons)<-TAM_pops
  # Convert to df
  DF_list[[ii]]<-plyr::ldply(ns.comparisons, cbind)
  #Addign column to identify parameter
  DF_list[[ii]]$Parameter<-rep(Stage_PAM[ii],nrow(DF_list[[ii]]))
}
# Summary of n.s. comparison with the respective population (single.modules[ii])
KW_results_DF<-do.call("rbind",DF_list)
print(KW_results_DF)
#