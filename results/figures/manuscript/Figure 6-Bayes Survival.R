# =============================================================================
# Name: Figure 5 
# Author: elopez
# Date: 29-08-2024
# Description: 
# TODO: 
# ============================================================================
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
#Abundance estimated by Bayes Prism
BP.theta <- readRDS(file = "D:/R_ScriptsPaper/Def_Objects/BP_Wu_Poles.rds")

# Population selection/order
TAM_pops <- c("FCN1","ISG15","CXCL9","IL1B","FOSB","FOLR2","APOE")
# Determining populations of interest
TAM_stringpattern <- grep(paste(TAM_pops,collapse = "|"), 
                          colnames(BP.theta), value = TRUE)

# % format
BP.theta <- BP.theta*100

#Subsetting TAM results (% scale)
TAM_estimation <- BP.theta[,TAM_stringpattern]

#Subsetting TAM with log10 scale of 100%
TAM_log_scale <- log10(BP.theta[,TAM_stringpattern])

# CLEANING clinical dataframe
# Import clinical data
clinical_ready <- readRDS(file = "D:/R_Scripts/EssentialObjects/clinical_ready.rds")
# Choosing clinical parameter to analyze
All_Param <-c("PAM50_Subtype","Pathologic_stage","PR_Status","ER_Status","Her2_Status")

# Check variables to clean
apply(clinical_ready[,All_Param], 2, table)

# Cleaning columns from clinical data 
to_clean<-c("NA","Normal","[Discrepancy]","[Not Available]","[Not Evaluated]","Indeterminate","Equivocal")

for (i in seq_along(to_clean)) {
  clinical_ready[All_Param][clinical_ready[All_Param]==to_clean[i]]<-NA
}
#
# ============================================================================
# Figure 5A: CELL TYPE CORRELATION MATRIX
# ============================================================================
#CORRELATION MATRIX
# All population (Only those with enough abundance)
library(Hmisc)
#res<-rcorr(BP.theta)
# Identifying  cell types with low abundances
# Calculate the xth percentile for each column in BP.theta
top_90th <- apply(BP.theta, 2, function(x) quantile(x, 0.995, na.rm = TRUE))

# Sort the variables based on the 75th percentile in decreasing order
sorted_top_90th <- names(top_90th)[order(top_90th, decreasing = TRUE)]

# 6 least abundant cell types
lowtypes <-tail(sorted_top_90th, 6)

res <- rcorr(BP.theta[,!colnames(BP.theta)%in%lowtypes])

# Perform hierarchical clustering
dist_matrix <- as.dist(1 - res$r)  # Convert correlation to distance
hclust_res <- hclust(dist_matrix, method = "complete")  # Hierarchical clustering

# Reorder the correlation matrix based on the clustering
ordered_r <- res$r[hclust_res$order, hclust_res$order]

# Convert the ordered correlation matrix to long format
cor_data <- reshape2::melt(ordered_r)

# Plot using ggplot2 with hierarchical clustering
celltype_corr_matrix <- ggplot(cor_data, aes(Var1, Var2, fill = value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(
    legend.title = element_text(size = 8),
    axis.text.x = element_text(angle = 45, size = 6,
                               vjust = 1, hjust = 1),
    axis.text.y = element_text(hjust = 1, size = 6),
    axis.title = element_blank()) +
  coord_fixed()
#
# ============================================================================
# SURVIVAL: DETERMINING CELLS WITH HIGH ESTIMATED ABUNDANCE
# ============================================================================
# Only abundant TAMs
TAMs_surv <- c("TAM_CXCL9","TAM_ISG15","TAM_FOLR2","TAM_APOE")
#-------REMOVE NA ROWS BEFORE PERFORMING INTERSECTION
patient_surv <- clinical_ready$patient[!is.na(clinical_ready$Evento)]

df_to_clean<-as.data.frame(TAM_estimation[rownames(TAM_estimation)%in%patient_surv,TAMs_surv])

#CHOOSE APPROPIATE LIMIT FOR EACH POPULATION OF INTEREST
apply(df_to_clean, 2, quantile,c(.75,.80,.90,.95,1))
#Select cells with % estimated abundance >=3.5
#Select in two fractions for better characterization
#Select cells with ES higher than 0.x for each column
Highest<-apply(df_to_clean,2,function(col) rownames(df_to_clean)[col>=4])
names(Highest)<-TAMs_surv

Snd.High<-apply(df_to_clean,2,function(col) rownames(df_to_clean)[col>=2 & col <4])
names(Snd.High)<-TAMs_surv

# For APOE TAMs we are lowering the threshold to 0.5
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
# Figure 5B: SURVIVAL-CXCL9 vs C1QA
# ============================================================================
#Defining abundant TAMs of interest
TAMs_surv <- c("TAM_CXCL9","TAM_ISG15","TAM_FOLR2","TAM_APOE")
# Selecting patients with high abundance of the phenotype of interest
my_list <- Final.list[TAMs_surv[!TAMs_surv%in%"TAM_ISG15"]]

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
curves_DF <- clinical_merge_abundance[,c("Months","Evento","Abundance")]

# Removing NA from 
curves_DF <- curves_DF[!is.na(curves_DF$Evento)&!is.na(curves_DF$Abundance),]

curves_DF$Abundance <- factor(curves_DF$Abundance, levels = c("TAM_CXCL9", 
                                                      "TAM_APOE", 
                                                      "TAM_FOLR2"))
#Survival analisis
a.km <- survfit(Surv(as.numeric(Months), Evento) ~ Abundance, data = curves_DF, type = "kaplan-meier")
#Para obtener median survival
#surv_median(a.km)
TAM_CXCL9_Curv <- ggsurvplot(fit = a.km, data = curves_DF, 
                             conf.int = FALSE, pval = F, 
                             surv.median.line = "none",
                             break.x.by = 50,
                             risk.table = "nrisk_cumevents",
                             risk.table.y.text.col = TRUE, 
                             risk.table.y.text = FALSE,
                             # Adjust this to make the numbers smaller
                             risk.table.fontsize = 2,  
                             legend.title = element_blank(),
                             legend.labs = c("High in CXCL9 TAMs",
                                             "High in APOE TAMs", 
                                             "High in FOLR2 TAMs"),
                             # Position the legend to the right
                             legend = "top",  
                             # Use theme_classic() to remove the grid
                             ggtheme = theme_classic() +           
                               theme(panel.grid = element_blank(),
                                     legend.text = element_text(size = 5),
                                     axis.title = element_text(size = 8),
                                     axis.text = element_text(size = 5),
                                     # Keep x and y axis lines
                                     axis.line = element_line(color = "black"))
)


# Perform pairwise log-rank tests to get pairwise p-values
pairwise_pvals <- pairwise_survdiff(Surv(as.numeric(Months), Evento) ~ Abundance, data = curves_DF)

# Extract the p-value matrix
print(pairwise_pvals$p.value)
#
# ============================================================================
# JOINING
# ============================================================================
# Extract the main plot and the risk table
Surv_plot <- TAM_CXCL9_Curv$plot
Surv_table <- TAM_CXCL9_Curv$table +
  theme(plot.title = element_text(size = 6))

# Combine the survival plot and the risk table into one plot
combined_surv_plot <- cowplot::plot_grid(Surv_plot, Surv_table, 
                                         ncol = 1, rel_heights = c(2, 0.8))

jpeg("D:/R_ScriptsPaper/Figures/Manuscrito/pre_Figure6.png", 
     res = 600,width = 18, height = 9, units = "cm")

# Combine the combined survival plot with the correlation matrix
cowplot::plot_grid(celltype_corr_matrix, combined_surv_plot,
                   ncol = 2, 
                   nrow = 1,
                   rel_widths = c(0.7, 1))
dev.off()
#