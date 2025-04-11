# ============================================================================
# Name: DEA FOSB vs ALL MONOCYTES
# Author: elopez
# Date: Ago-04-2024
# Description: 
# TODO: 
# ============================================================================
library(Seurat)
library(ggplot2)

# Object with phenotypes assigned in all tissue types
Tissue_Pheno <-readRDS(file = "D:/R_ScriptsPaper/Def_Objects/Tissue_Pheno.rds")

# Only Periphery
BLOOD <- subset(Tissue_Pheno, subset= Tissue == "Blood")

# Identifying high FOSB cells in blood (based on median score value)
blood_brca_cells <- BLOOD@meta.data %>%
  filter(TissueType == "Blood: BRCA") %>%
  # Compute the 75th percentile of FOSB_score within this subset
  mutate(quantile_75 = quantile(FOSB_score, 0.5, na.rm = TRUE)) %>%
  # Select rows where FOSB_score is greater than the 75th percentile
  filter(FOSB_score > quantile_75) %>%
  # Select the rownames of the filtered rows
  rownames()

# Rename identity of those cells
BLOOD <- SetIdent(object = BLOOD, cells = blood_brca_cells, value = "Mono_FOSB")
# Second group, all of the other monocytes
other_cells <- setdiff(colnames(BLOOD), blood_brca_cells)
BLOOD <- SetIdent(object = BLOOD, cells = other_cells, value = "Others")
# Differential expression
BloodDF <- FindMarkers(BLOOD, logfc.threshold = 0.1, min.pct = 0.05,
                     #Monocitos FOSB en BRCA: Blood (+)
                     ident.1 = "Mono_FOSB",
                     #Resto de poblaciones de monos (-)
                     ident.2 = "Others",
                     #Datos normalizados
                     assay = "RNA", slot = "data")

library(dplyr)
# Define patterns to remove structural genes
bad_patterns <- 'RP[0-9]|^MT|[0-9]{4}|^RPL|^RPS|[0-9]orf[0-9]|^ATP'

filtered_bloodDF <- BloodDF %>%
  # Filter out rows with bad patterns in rownames and significant p-values
  filter(!grepl(bad_patterns, rownames(BloodDF)), p_val_adj < 0.01) %>%
  # Sort by avg_log2FC in descending order
  arrange(desc(avg_log2FC))


# Save results with all genes
write.csv(filtered_bloodDF, file = "D:/R_ScriptsPaper/Def_Objects/DEG_FOSBvsAll.csv", 
          row.names = TRUE)
