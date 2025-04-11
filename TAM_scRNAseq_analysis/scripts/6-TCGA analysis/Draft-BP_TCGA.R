# =============================================================================
# Name: BAYESPRISM-ANALYSIS OF ESTIMATED ABUNDANCE 
# Author: elopez
# Date: 06-08-2024
# Description: 
# TODO: 
# ============================================================================
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
# CLEANING clinical dataframe
# Import clinical data
clinical_ready <- readRDS(file = "D:/R_Scripts/EssentialObjects/clinical_ready.rds")
# Choosing clinical parameter to analyze
All_Param<-c("PAM50_Subtype","Pathologic_stage","PR_Status","ER_Status","Her2_Status")

# Check variables to clean
apply(clinical_ready[,All_Param], 2, table)

# Cleaning columns from clinical data 
to_clean<-c("NA","Normal","[Discrepancy]","[Not Available]","[Not Evaluated]","Indeterminate","Equivocal")

for (i in seq_along(to_clean)) {
  clinical_ready[All_Param][clinical_ready[All_Param]==to_clean[i]]<-NA
}

#Abundance estimated by Bayes Prism
#BP.theta <- readRDS(file = "D:/R_ScriptsPaper/Objects/Temporal/BayesPrism/theta_MainPheno_HSP.rds")
BP.theta <- readRDS(file = "D:/R_ScriptsPaper/Objects/Temporal/BayesPrism/theta_Polos.rds")

# Population selection/order
TAM_pops<-c("FCN1","ISG15","CXCL9","FOSB","HSPs","HLA_II","C1QA","APOE")

# % format
BP.theta<-BP.theta*100

# Determining populations of interest
TAM_stringpattern <- grep(paste(TAM_pops,collapse = "|"), 
                colnames(BP.theta), value = TRUE)

#Subsetting TAM results (% scale)
TAM_estimation<-BP.theta[,TAM_stringpattern]

#Subsetting TAM with log10 scale of 100%
TAM_log_scale<-log10(BP.theta[,TAM_stringpattern])
#
# ============================================================================
# Section 1: EVALUATING TAM ABUNDANCE
# ============================================================================
# Visualize co existance within patient
library(pheatmap)
library(gplots)
#Palette color
spectro<-RColorBrewer::brewer.pal(11, "Spectral")
# Aesthetic
TAM_est_DF<-TAM_estimation
colnames(TAM_est_DF)<-paste0("TAM_",colnames(TAM_est_DF))
# Heatmap of %
pheatmap(t(TAM_est_DF), col=rev(spectro), cluster_cols=T, cluster_rows=T, 
         fontsize_col=6, fontsize_row = 8,
         border_color=NA, 
         show_colnames = F,
         annotation_legend = F, 
         annotation_names_row = F)

# ABUNDANCE (%) OF TAM POPS
Total_TAMs<-as.data.frame(TAM_estimation)
Total_TAMs$Total_TAM<-apply(Total_TAMs,1,sum)
# Median value
median_value <- apply(Total_TAMs, 2, median)
# Setting order by median_value
orden_abun<-names(median_value)[order(median_value, decreasing = T)]
# Format for plotting
gg_Total_TAMs<-reshape2::melt(Total_TAMs)

# Custom in X axis
original_labels <- unique(gg_Total_TAMs$variable) 
new_labels <- paste0("TAM_",original_labels) 
# Create a named vector for new labels 
xLabels <- setNames(new_labels, original_labels) 

ggplot(gg_Total_TAMs, aes(x = variable, y = value)) +
  geom_boxplot()+
  #Orden eje x
  theme_minimal()+
  scale_x_discrete(limits=orden_abun,
                   labels=xLabels)+
  ylab("Percentage within the TME")+
  theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1),
        axis.title.x = element_blank())

#CORRELATION FOR ONLY TAM
#PerformanceAnalytics::chart.Correlation(TAM_estimation, histogram=TRUE, pch=19)
#
# ============================================================================
# Section 2: EXPLORING ABUNDANCES IN THE TME
# ============================================================================
#Median values
medians <- apply(BP.theta, 2, median)
sorted_medinas<-names(medians)[order(medians, decreasing = T)]
# Reshape the dataframe to long format
df_long <- reshape2::melt(BP.theta)

# Percentages in the TME
ggplot(df_long, aes(x = Var2, y = value)) +
  geom_boxplot()+
  #Orden eje x
  theme_minimal()+
  scale_x_discrete(limits=sorted_medinas)+
  labs(title="Population abundance in TME")+
  theme(axis.text.x = element_text(angle = 45,hjust=1, vjust=1),
        axis.title = element_blank())

#CORRELATION MATRIX
# All population (Only those with enough abundance)
library(Hmisc)
res<-rcorr(BP.theta)

# Perform hierarchical clustering
dist_matrix <- as.dist(1 - res$r)  # Convert correlation to distance
hclust_res <- hclust(dist_matrix, method = "complete")  # Hierarchical clustering

# Reorder the correlation matrix based on the clustering
ordered_r <- res$r[hclust_res$order, hclust_res$order]

# Convert the ordered correlation matrix to long format
cor_data <- reshape2::melt(ordered_r)

# Plot using ggplot2 with hierarchical clustering
ggplot(cor_data, aes(Var1, Var2, fill = value)) + 
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(hjust = 1),
        axis.title = element_blank()) +
  coord_fixed()
#
# ============================================================================
# Section 3.1: CLINICAL ASSOCIATIONS-STATISTICAL TESTS
# ============================================================================
# Merge the data frames, using row names of TAM_estimation as the patient column
clin_TAMestim <- merge(clinical_ready[,c("patient",All_Param,"Evento")], TAM_estimation, 
                       by.x = "patient", by.y = "row.names")
#-------OPTION: ONLY ROWS WITH SURVIVAL DATA
# Filter out rows with NA in Evento
#clin_TAMestim <- clin_TAMestim %>%
 # filter(!is.na(Evento))
#-------END

# ----------SHAPIRO TEST (NORMALITY TEST) FOR TAM ABUNDANCE
for (hh in seq_along(All_Param)) {
  # Subsetting df for analysis
  sub_df<-clin_TAMestim[,c(All_Param[hh],TAM_pops)]
  colnames(sub_df)<-c("Parameter",TAM_pops)
  # spliting by clinical parameter and remove rows with NA
  split_df <- sub_df %>%
    filter(!is.na(Parameter)) %>%
    group_split(Parameter)
  
  # List to store shapiro results
  result_df<-list()
  # Testing distribution of each parameter
  for (ii in seq_along(split_df)) {
    # Selecting respective df of each parameter (eliminating parameter column)
    shap.df<-split_df[[ii]][-1]
    #Shapiro test for each signature per clinical parameter
    shapiro_results <- apply(shap.df, 2, shapiro.test)
    # Pval for each score
    shap_res<-vector()
    for (jj in seq_along(shapiro_results)) {
      # If pval<0.05, data is not normally distributed
      shap_res[jj]<-shapiro_results[[jj]]$p.value>0.05
    }
    # A result_df for each parameter = length(split_df)
    result_df[[ii]]<-data.frame(names(shapiro_results),shap_res)
    colnames(result_df[[ii]])<-c(unique(split_df[[ii]]$Parameter), "Normal_distribution")
  }
  print(All_Param[hh])
  print(result_df)
}
# ----------KRUSKAL-WALLIS for 3 or more groups; non-normal distribution
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

#------------MANN WHITNEY WILCOXON
# Define the columns of interest
Pos_Neg <- c("PR_Status", "ER_Status", "Her2_Status")
# Subset the dataframe
wlx_df <- clin_TAMestim[, c(Pos_Neg, TAM_pops)]
# Initialize an empty list to store results
mann_whitney_results <- list()
# Apply the Mann-Whitney U test for each combination of clinical status and TAM population
for (clinical_status in Pos_Neg) {
  for (tam_pop in TAM_pops) {
    test_result <- wilcox.test(wlx_df[[tam_pop]] ~ wlx_df[[clinical_status]], data = wlx_df)
    mann_whitney_results[[paste(clinical_status, tam_pop, sep = "_")]] <- test_result$p.value
  }
}
# Convert results list to a dataframe and transpose
MannWhit_df <- as.data.frame(t(as.data.frame(mann_whitney_results)))
colnames(MannWhit_df) <- c("P_Value")
MannWhit_df$Significant<-MannWhit_df$P_Value<0.01
# Print the results
print(MannWhit_df[MannWhit_df$Significant,])
#
# ============================================================================
# Section 3.2: CLINICAL ASSOCIATIONS-BOXPLOT
# ============================================================================
# Merge the data frames, using row names of TAM_log_scale as the patient column
clin_TAMlog <- merge(clinical_ready[,c("patient",All_Param,"Evento")], TAM_log_scale, 
                       by.x = "patient", by.y = "row.names")

#----------OPTION 2: ONLY ROWS WITH SURVIVAL DATA
# Filter out rows with NA in Evento and in PAM50_Subtype
#clin_TAMlog <-clin_TAMlog[!is.na(clin_TAMlog$Evento),]
#----------END

# Eliminating rows with no estimated abundance for TAM populations
long_clinTAM <- clin_TAMlog %>%
  pivot_longer(cols = HSPs:CXCL9,    # Specify the range of columns to pivot
               names_to = "variable",
               values_to = "value") 

# Setting order for plotting
# Adding TAM prefix
long_clinTAM$variable<-paste0("TAM_",long_clinTAM$variable)
# Setting TAM population order
long_clinTAM$variable<-factor(long_clinTAM$variable,
                              levels = paste0("TAM_",TAM_pops))
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
# Define the layout matrix with an empty space in the first row, third column
layout <- rbind(c(1, 2, NA),  # First row with two plots and one empty space
                c(3, 4, 5))   # Second row with three plots
# Arrange the plots according to the layout
gridExtra::grid.arrange(grobs = boxpl_list, layout_matrix = layout)
#
# ============================================================================
# Section 4: HAZARD RATIOS
# ============================================================================
# FCN1 AND HSPs ARE NOT ABUNDANT ENOUGH FOR THIS ANALYSIS. THEY WONT BE CONSIDERED
# Moving rownames before merging
patient_TAM_abundance<-tibble::rownames_to_column(as.data.frame(TAM_estimation), 
                                                  var = "patient")
# Merging with clinical data
clinical_TAM_abundance<-dplyr::left_join(clinical_ready,patient_TAM_abundance,by="patient")

# Selecting only necessary columns for survival analysis
Hazard_ratio<-clinical_TAM_abundance[,c("Months","Evento",TAM_stringpattern)]
# Removing rows with NA in Evento column
Hazard_ratio<-Hazard_ratio[!is.na(Hazard_ratio$Evento),]
#Cox model with desired populations
overall_TAM <- coxph(Surv(Months, Evento) ~ ISG15+ CXCL9+ FOSB+
                       HLA_II+ C1QA+ APOE, data = Hazard_ratio)

# Create a forest plot of the Cox model
survminer::ggforest(overall_TAM, data = Hazard_ratio)
#Resultados
summary(overall_TAM)$coef
#pval
cox.zph(overall_TAM)
# ============================================================================
# Section 5: SURVIVAL CURVES
# ============================================================================
#CHOOSE APPROPIATE LIMIT FOR EACH POPULATION OF INTEREST
apply(TAM_estimation, 2, quantile,c(.75,.80,.90,.95,1))

# Selecting patients with abundance >5% for population of interest 
# (ISG15 is not associated with prognosis. High HLAII patients with survival data is a small group)
TAMs_surv <-c("C1QA","CXCL9")
High_Patients<-list()

for (i in 1:length(TAMs_surv)) {
  if(TAMs_surv[i]%in%c("HLA_II","FOSB")){
    #0.1%
    limit<-0.5
  }else{
    #0.5%
    limit<-0.5
  }
  #First remove NA for the respective column
  cleaned_df <- clinical_TAM_abundance[!is.na(clinical_TAM_abundance[,TAMs_surv[i]]),]
  #Select patients above the limit
  boslog <- cleaned_df[,grep(TAMs_surv[i],colnames(cleaned_df))]>limit
  High_Patients[[i]] <- cleaned_df$patient[boslog]
}
names(High_Patients)<-TAMs_surv

library(venn)
HIGH_venn<-venn(High_Patients)
#Extrayendo interseccion con funcion attribute
HIGH_inter<-attr(x = HIGH_venn, "intersections")

# Selected combination
my_list<-HIGH_inter[names(HIGH_inter)%in%TAMs_surv]

# Second option. Top 6 phenotypes
#my_list<-HIGH_inter[names(sort(lengths(HIGH_inter), decreasing = T))[1:4]]

# Convert the list to a named vector
abundance_subtitution <- unlist(my_list)
# Adjust the names
names(abundance_subtitution) <- paste0("H_",rep(names(my_list), lengths(my_list)))

#name vector as dataframe
abund_df<-data.frame(patient= as.vector(abundance_subtitution),
                     Abundance= names(abundance_subtitution))
# Joining df
clinical_merge_abundance<-merge(clinical_TAM_abundance,abund_df,by='patient',all=TRUE)

# Only necessary columns
curves_DF<-clinical_merge_abundance[,c("Months","Evento","Abundance")]

# Removing NA from evento column
curves_DF<-curves_DF[!is.na(curves_DF$Evento)&!is.na(curves_DF$Abundance),]

#Survival analisis
a.km <- survfit(Surv(as.numeric(Months), Evento) ~ Abundance, data = curves_DF, type = "kaplan-meier")
#Para obtener median survival
#surv_median(a.km)
ggsurvplot(fit=a.km, data = curves_DF, 
           palette = c("red","blue"),
           conf.int=FALSE, pval=T, surv.median.line = "none",break.x.by=50,
           risk.table = "nrisk_cumevents",
           risk.table.y.text.col = T, # colour risk table text annotations.
           risk.table.y.text = FALSE,
           legend.title=element_blank(),xlab = NULL,ylab = NULL)
# Perform pairwise log-rank tests to get pairwise p-values
pairwise_pvals <- pairwise_survdiff(Surv(as.numeric(Months), Evento) ~ Abundance, data = curves_DF)

# Extract the p-value matrix
print(pairwise_pvals$p.value)
#
# ============================================================================
# Irrelevant: CXCL9 Stage I-II vs III-IV
# ============================================================================
# Merge the data frames, using row names of TAM_log_scale as the patient column
clin_TAMlog <- merge(clinical_ready[,c("patient",All_Param,"Evento")], TAM_log_scale, 
                     by.x = "patient", by.y = "row.names")

#CXCL9
stages_df<-clin_TAMlog[,c("Pathologic_stage","CXCL9")]
library(dplyr)
library(ggplot2)
library(ggpubr)

# Combine stages
stages_df <- stages_df %>%
  mutate(Pathologic_stage = case_when(
    Pathologic_stage %in% c("Stage_I", "Stage_II") ~ "Stage_I_II",
    Pathologic_stage %in% c("Stage_III", "Stage_IV") ~ "Stage_III_IV",
    TRUE ~ Pathologic_stage
  )) %>%
  filter(!is.na(Pathologic_stage))  # Remove rows with NA in Pathologic_stage

# Remove rows below log10(1) = 0 in the CXCL9 column
stages_df_filtered <- stages_df %>%
  filter(CXCL9 >= 0)

# Revert CXCL9 to the original scale for the Mann-Whitney test
stages_df_filtered$CXCL9_original <- 10^stages_df_filtered$CXCL9

# Perform the Mann-Whitney test between the two groups
mann_whitney_test <- wilcox.test(CXCL9_original ~ Pathologic_stage, data = stages_df_filtered)
p_value <- mann_whitney_test$p.value

# Create the boxplot
ggplot(stages_df_filtered, aes(x = Pathologic_stage, y = CXCL9, fill = Pathologic_stage)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Pathologic Stage", y = "CXCL9 (log10 scale)") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.y = max(stages_df_filtered$CXCL9) + 0.1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle(paste("Mann-Whitney p-value:", signif(p_value, digits = 3)))

#
# ============================================================================
# MYELOID LYMPHOID 
# ============================================================================
# Define the myeloid and lymphoid cell types
myeloid_cells <- c("FCN1", "ISG15", "CXCL9", "FOSB", "HSPs", "HLA_II", 
                   "C1QA", "APOE", "cDC1", "cDC2", "Mature DC", "pDC")

lymphoid_cells <- c("B cell", "Plasmatic B cell", "T CD4", "T CD8", 
                    "NK cell", "T regs")

# Calculate the row sums for myeloid cells
myeloid_sums <- rowSums(BP.theta[, myeloid_cells], na.rm = TRUE)

# Calculate the row sums for lymphoid cells
lymphoid_sums <- rowSums(BP.theta[, lymphoid_cells], na.rm = TRUE)

# Create a new dataframe with the results of the sums
TILs_sums <- data.frame(
  Myeloid_Sum = myeloid_sums,
  Lymphoid_Sum = lymphoid_sums
)
#
#%%%: CORRELACION de ES con PERCENT####
TCGA_ES<-readRDS(file = "D:/R_Scripts/EssentialObjects/TCGA_enrich.rds")
colnames(TCGA_ES)[3]<-"CXCR4"
colnames(TCGA_ES)<-paste(colnames(TCGA_ES),"_ES",sep = "")
TCGA_ES<-tibble::rownames_to_column(as.data.frame(TCGA_ES), var = "patient")

TAM_estimation<-BP.theta[,grep("TAM_",colnames(BP.theta))]
TAM_estimation<-tibble::rownames_to_column(as.data.frame(TAM_estimation), var = "patient")
#Juntando con porcentajes
ES_Por<-merge(TCGA_ES,TAM_estimation,by='patient',all=TRUE)

Pobs<-c("HLA_II","C1QA","SPP1","ISG15","FCN1","CXCR4")
lista<-list()
for (i in 1:length(Pobs)) {
  ESP<-ES_Por[,grep(Pobs[i],colnames(ES_Por))]
  colnames(ESP)<-c("ES","Percentage")
  lista[[i]]<-ggplot(ESP, aes(x=ES, y=Percentage)) + 
    geom_point()+
    geom_smooth()+
    labs(title = paste("TAM_",Pobs[i],sep = ""))+
    xlab("GSVA: Enrichment score")+
    ylab("BayesPrism: Percentage")+
    theme(plot.title = element_text(hjust = 0.5))
}
ggpubr::ggarrange(plotlist=lista,ncol = 3, nrow = 2)
#
#%%%: TREND ANALYSIS####
# Install and load the ordinal package
#install.packages("ordinal")
#Prueba con datos sin logaritmo
TAM_nolog<-BP.theta[,grep("TAM_",colnames(BP.theta))]
TAM_nolog<-tibble::rownames_to_column(as.data.frame(TAM_nolog), var = "patient")
px_TAM_logscale<-reshape2::melt(TAM_nolog)

df<-merge(clinical_ready,px_TAM_logscale,by='patient',all=TRUE)

A<-df[,c("variable","value","Pathologic_stage")]
colnames(A)<-c("Module","ES","Parameter")
dfdf<-A[!A$Parameter%in%c("NA/Other","Stage X","NA"),]
#Eliminando na 
dfdf<-dfdf[!is.na(dfdf$ES),]
#Elegir Poblacion
Tend_ready<-dfdf[dfdf$Module=="TAM_HLA_II",]
Tend_ready$Parameter<-factor(Tend_ready$Parameter, 
                       levels = c("Stage_I","Stage_II","Stage_III","Stage_IV"))

library(ordinal)
# Fit proportional odds ordinal logistic regression
model <- clm(Parameter ~ ES, data = Tend_ready)

# Summary of the model
summary(model)

#
#GRAIFCOS DE COMBINACIONES----------
#HEATMAP
TCGA_enrich<-readRDS(file="D:/R_Scripts/EssentialObjects/TCGA_enrich.rds")
res<-apply(TCGA_enrich, 2,as.numeric)
colnames(res)<-paste(colnames(res),"score", sep = " ")
library(pheatmap)   
library(gplots)
#Palette color
Rojos<-RColorBrewer::brewer.pal(9, "YlOrRd")
#Definiendo colores de signatures
F.DF<-as.data.frame(letters[1:6])
rownames(F.DF)<-colnames(res)
colnames(F.DF)<-"Signatures"
ann_colors <- list(Signatures = 
                     c(a = "brown", b = "green", c = "turquoise",
                       d = "blue",e = "#FFCC00",f = "black"))
C<-pheatmap(t(res), col=Rojos, cluster_cols=T, cluster_rows=T, 
            fontsize_col=6, fontsize_row = 8, legend_labels = "NES",
            border_color=NA, angle_col = 90, annotation_row = F.DF,
            annotation_legend = F, annotation_names_row = F,
            annotation_colors = ann_colors)
#Convertir a ggplot para usarlos con patchwork
HEAT_TCGA<-ggplotify::as.ggplot(C)
#

#HEATMAP_HIGH
ORDEN<-c("HLA_II","C1QA","SPP1","FCN1","VEGFA","ISG15")
LowTCGA<-readRDS(file = "D:/R_Scripts/LowTCGA.rds")
DF_CR<-clinical_ready[,c("patient","TAM_Pob")]
DF_CRLogic<-DF_CR
#Indicando si paciente es high para alguna poblacion
for (i in 1:length(ORDEN)) {
  A<-DF_CR[grepl(ORDEN[i],DF_CR$TAM_Pob),c("patient")]
  #Anadiendo valor a 
  DF_CR[ORDEN[i]]<-ifelse(DF_CR$patient%in%A,"High",
                          ifelse(DF_CR$patient%in%LowTCGA[[ORDEN[i]]],"Low",
                                 "Medium"))
  DF_CRLogic[ORDEN[i]]<-ifelse(DF_CR$patient%in%A,6,
                               ifelse(DF_CR$patient%in%LowTCGA[[ORDEN[i]]],1,
                                      0))
}
WTFALL<-DF_CR[,c("patient",ORDEN)]
#Para fijar orden.
WTFALL$suma <- apply(DF_CRLogic[,-c(1:2)], 1, sum)
WTFALL<-WTFALL[order(WTFALL$suma, decreasing = T),]
WTFALL$suma<-NULL
DF<-reshape2::melt(WTFALL, id.vars = "patient")
#Ordenado por high, low o medium
DF<-DF[order(DF$value),]
#Orden de pacientes
DF$patient<-factor(DF$patient, levels = unique(DF$patient))
#Orden de poblaciones
DF$variable<-factor(DF$variable, levels = rev(ORDEN))
HEAT_Discrete<-ggplot(DF,aes(patient,variable,fill=as.factor(value)))+
  geom_tile()+
  #Resetea el fondo del grafico
  theme_classic()+
  #Colores de etiquetas
  scale_fill_manual(values = c("#009900","purple","grey"))+
  ylab(NULL)+xlab("BRCA patients (n=1099)")+
  #Nombre de ejes y etiquetas
  labs(#title = "Patients high in TAMs",
       fill = "TAM enrichment")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold"),
        #Estilo de los paneles
        panel.border = element_rect(fill=NA),
        plot.title = element_text(face = "bold",size = 12,hjust = 0.5, family = "serif", colour = "blue"),
        #La aplitud del eje X es ? veces mayor que Y y/x
        aspect.ratio = 0.3)
#
#PASTEL
#Preparando datos
CAKKE<-stringr::str_count(clinical_ready$TAM_Pob,":")
CAKKE<-ifelse(is.na(CAKKE),"Mid/Low",CAKKE)

CAKKE<-replace(CAKKE,CAKKE%in%c("3","4","5"),"4 o more TAM population")

DF_CAKKE<-as.data.frame(table(CAKKE)/length(CAKKE)*100)
DF_CAKKE["dummy"]<-"dummy"
#Cambiando a mano
DF_CAKKE$CAKKE<-c(paste(c("1","2","3","4 o more"),"TAM populations",sep = " "),"None")

library(scales)
library(ggrepel)
#BACKGROUND BLANCO
blank_theme <- theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        #Sin texto en los tick marks
        axis.text = element_blank(),
        panel.border = element_blank(),
        legend.title =element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 8),
        panel.grid=element_blank())

Pastel_TAM<-ggplot(DF_CAKKE, aes(x=dummy,y=Freq, fill=CAKKE))+
  geom_col(colour="black")+
  blank_theme +#labs(title = "")+
  #Colores de las fracciones
  labs(fill="Patients enriched with")+
  #El parametro x=1.55 dentro de aes(), posiciona el numero fuera del pie
  ggrepel::geom_label_repel(aes(x = 1.5,label = percent(Freq/100, accuracy = 0.1)),
                            position = position_stack(vjust = 0.5), size=2.5, 
                            #Quita texto de los cuadros de colores de la leyenda
                            show.legend = F)+
  theme(plot.title = element_text(hjust = 0.5,
                                  #colour=paleta[i], 
                                  size = rel(0.8)),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = rel(0.6)))+
  coord_polar("y")
#
#FREQUENTLY FOUND TOGETHER
#BULK Pairs
G<-stringr::str_count(clinical_ready$TAM_Pob,":")
Pobs<-c("FCN1","C1QA","VEGFA","ISG15","HLA_II","SPP1")
#Todas las combinaciones posibles entre 2
pb.df<-t(combn(Pobs, 2))
#Las combinaciones con 2 o mas poblaciones
CombFrec<-clinical_ready$TAM_Pob[G>=1]
#Contando las combinaciones
Rloop<-list()
Rloop.name<-vector()
for (i in 1:nrow(pb.df)){
  #La combinacion a buscar
  patt<-paste("(?=.*",pb.df[i,1],")(?=.*",pb.df[i,2],")",sep = "")
  #Cuantos de cada combinacion
  Rloop[[i]]<-length(CombFrec[grepl(patt,CombFrec, perl = T)])
  Rloop.name[i]<-paste(pb.df[i,1],pb.df[i,2],sep = ":")
}
names(Rloop)<-Rloop.name
#De mayor a menor
Rloop<-Rloop[order(unlist(Rloop), decreasing = T)]
Rl<-plyr::ldply(Rloop, rbind)
colnames(Rl)<-c("Pairs","Num.Patients")
Rl$Pairs<-factor(Rl$Pairs, levels = rev(unique(Rl$Pairs)))
BULK_Pairs<-ggplot(Rl, aes(x=Num.Patients, y=Pairs))+
  geom_point(colour="red",size=2)+
  #geom_segment(aes(xend=CC), yend=0) +
  expand_limits(x=150) +
  ylab("Co-expressing signatures")+
  xlab("Number of BRCA patients")+
  labs(title = "Co-expressed signature pairs")+
  geom_text(label=Rl$Num.Patients,hjust=-1,size=2.3)+
  theme_classic()+
  theme(axis.title = element_text(size = rel(0.8)),
        axis.text.y=element_text(size = rel(0.7)),
        panel.border = element_rect(fill=NA))
#
