library (SingleR)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(gridExtra)
library(reshape2)

#load seurat objects
setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/")
load("myel.RData")
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/")
load("macro_7dpi.RData")
macro_B <- macro
rm(macro)
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
load("macro_brennan.RData")
macro_C <- macro
rm(macro)

##AxB -----------------------------------------------------------------
#prepare reference labels
reference_labels <- myel@active.ident

#prepare query clusters
query_clusters <- as.factor(macro_B@active.ident)

#Run Single R and extract results
singleR_result_AxB <- SingleR(test = GetAssayData(macro_B, slot = "data"), 
                          ref = GetAssayData(myel, slot = "data"), 
                          labels = reference_labels)
predicted_cell_types <- singleR_result$labels
macro_B$predicted_cell_types <- predicted_cell_types

##Compare how clusters in the query correspond to predicted cell types
# Create a contingency table of query clusters vs predicted cell types
contingency_table <- table(query_clusters = query_clusters, 
                           predicted_cell_types = macro_B$predicted_cell_types)
df_long <- melt(contingency_table)

ggplot(df_long, aes (x = "query_clusters", y = "predicted_cell_types", fill = "Freq")) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(x = "Dataset A Cluster", y = "Dataset B Cluster", title = "Dataset A x Dataset B", fill = "Frequency of predicted cells from Dataset B") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




pheatmap(contingency_table, 
         main = "Dataset A vs Dataset B", 
         color = colorRampPalette(c("white", "blue"))(50), 
         scale = "none",
         cluster_cols = FALSE,
         cluster_rows = FALSE) 
##AxC -----------------------------------------------------------------
#load seurat objects
setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/")
load("myel.RData")
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
load("macro_brennan.RData")

#prepare reference labels
reference_labels <- myel@active.ident

#prepare query clusters
query_clusters <- as.factor(macro@active.ident)

#Run Single R and extract results
singleR_result <- SingleR(test = GetAssayData(macro, slot = "data"), 
                          ref = GetAssayData(myel, slot = "data"), 
                          labels = reference_labels)
predicted_cell_types <- singleR_result$labels
macro$predicted_cell_types <- predicted_cell_types

##Compare how clusters in the query correspond to predicted cell types
# Create a contingency table of query clusters vs predicted cell types
contingency_table <- table(query_clusters = query_clusters, 
                           predicted_cell_types = macro$predicted_cell_types)

pheatmap(contingency_table, 
         main = "Dataset A vs Dataset C", 
         color = colorRampPalette(c("white", "blue"))(50), 
         scale = "none",
         cluster_cols = FALSE,
         cluster_rows = FALSE) 
##BxC -----------------------------------------------------------------
#load seurat objects
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/")
load("macro_7dpi.RData")
macro_B <- macro
rm(macro)
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
load("macro_brennan.RData")
macro_C <- macro
rm(macro)

#prepare reference labels
reference_labels <- macro_B@active.ident

#prepare query clusters
query_clusters <- as.factor(macro_C@active.ident)

#Run Single R and extract results
singleR_result <- SingleR(test = GetAssayData(macro_C, slot = "data"), 
                          ref = GetAssayData(macro_B, slot = "data"), 
                          labels = reference_labels)
predicted_cell_types <- singleR_result$labels
macro_C$predicted_cell_types <- predicted_cell_types

##Compare how clusters in the query correspond to predicted cell types
# Create a contingency table of query clusters vs predicted cell types
contingency_table <- table(query_clusters = query_clusters, 
                           predicted_cell_types = macro_C$predicted_cell_types)

pheatmap(contingency_table, 
         main = "Dataset B vs Dataset C", 
         color = colorRampPalette(c("white", "blue"))(50), 
         scale = "none",
         cluster_cols = FALSE,
         cluster_rows = FALSE) 
