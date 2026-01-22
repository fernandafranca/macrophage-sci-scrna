library (SingleR)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(reshape2)
library(patchwork)
library(tidyr)
library(dplyr)

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
reference_labels_AxB <- myel@active.ident

#prepare query clusters
query_clusters_AxB <- as.factor(macro_B@active.ident)

#Run Single R and extract results
singleR_result_AxB <- SingleR(test = GetAssayData(macro_B, slot = "data"), 
                          ref = GetAssayData(myel, slot = "data"), 
                          labels = reference_labels_AxB,
                          de.method = "wilcox")
predicted_cell_types_AxB <- singleR_result_AxB$labels
macro_B$predicted_cell_types <- predicted_cell_types_AxB


##Compare how clusters in the query correspond to predicted cell types
# Create a contingency table of query clusters vs predicted cell types
contingency_table_AxB <- table(query_clusters = query_clusters_AxB, 
                           predicted_cell_types = macro_B$predicted_cell_types)
rownames(contingency_table_AxB) <- paste0("Cluster ", rownames(contingency_table_AxB))


# Convert to data frame, ensure row names are stored properly
df_long_AxB <- contingency_table_AxB %>%
  as.data.frame() %>% 
  rename(value = "Freq") %>%  # Rename "Freq" column to "value"
  mutate(across(where(is.factor), as.character))  # Convert any factors to characters

heatmap1 <- ggplot(df_long_AxB, aes(x = predicted_cell_types, y = query_clusters, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#F4A3C0") +
  labs(title = "Macrophage Cluster Matching (Dataset A vs Dataset B)",
       x = "Predicted Cell Types (Dataset A)", 
       y = "Query Clusters (Dataset B)",
       fill = "Number of cells predicted on test dataset") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
heatmap1

setwd("~/Google Drive/My Drive/Gensel lab/R studies/SingleR/")
save(singleR_result_AxB, file = "singleR_AxB_wilcox.RData")
save(df_long_AxB, file = "df_long_AxB_wilcox.RData")

##AxC -----------------------------------------------------------------
#prepare reference labels
reference_labels_AxC <- myel@active.ident

#prepare query clusters
query_clusters_AxC <- as.factor(macro_C@active.ident)

#Run Single R and extract results
singleR_result_AxC <- SingleR(test = GetAssayData(macro_C, slot = "data"), 
                          ref = GetAssayData(myel, slot = "data"), 
                          labels = reference_labels_AxC,
                          de.method = "wilcox")
predicted_cell_types_AxC <- singleR_result_AxC$labels
macro_C$predicted_cell_types <- predicted_cell_types_AxC

##Compare how clusters in the query correspond to predicted cell types
# Create a contingency table of query clusters vs predicted cell types
contingency_table_AxC <- table(query_clusters = query_clusters_AxC, 
                           predicted_cell_types = macro_C$predicted_cell_types)
rownames(contingency_table_AxC) <- paste0("Cluster ", rownames(contingency_table_AxC))

# Convert to data frame, ensure row names are stored properly
df_long_AxC <- contingency_table_AxC %>%
  as.data.frame() %>% 
  rename(value = Freq) %>%  # Rename "Freq" column to "value"
  mutate(across(where(is.factor), as.character))  # Convert any factors to characters

heatmap2 <- ggplot(df_long_AxC, aes(x = predicted_cell_types, y = query_clusters, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#F4A3C0") +
  labs(title = "Macrophage Cluster Matching (Dataset A vs Dataset C)",
       x = "Predicted Cell Types (Dataset A)", 
       y = "Query Clusters (Dataset C)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
heatmap2

setwd("~/Google Drive/My Drive/Gensel lab/R studies/SingleR/")
save(singleR_result_AxC, file = "singleR_AxC_wilcox.RData")
save(df_long_AxC, file = "df_long_AxC_wilcox.RData")


##BxC -----------------------------------------------------------------
#prepare reference labels
reference_labels_BxC <- macro_B@active.ident

#prepare query clusters
query_clusters_BxC <- as.factor(macro_C@active.ident)

#Run Single R and extract results
singleR_result_BxC <- SingleR(test = GetAssayData(macro_C, slot = "data"), 
                          ref = GetAssayData(macro_B, slot = "data"), 
                          labels = reference_labels_BxC,
                          de.method = "wilcox")
predicted_cell_types_BxC <- singleR_result_BxC$labels
macro_C$predicted_cell_types <- predicted_cell_types_BxC

##Compare how clusters in the query correspond to predicted cell types
# Create a contingency table of query clusters vs predicted cell types
contingency_table_BxC <- table(query_clusters = query_clusters_BxC, 
                           predicted_cell_types = macro_C$predicted_cell_types)

rownames(contingency_table_BxC) <- paste0("Cluster ", rownames(contingency_table_BxC))

# Convert to data frame, ensure row names are stored properly
df_long_BxC <- contingency_table_BxC %>%
  as.data.frame() %>% 
  rename(value = Freq) %>%  # Rename "Freq" column to "value"
  mutate(across(where(is.factor), as.character))  # Convert any factors to characters

heatmap3 <- ggplot(df_long_BxC, aes(x = predicted_cell_types, y = query_clusters, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#F4A3C0") +
  labs(title = "Macrophage Cluster Matching (Dataset B vs Dataset C)",
       x = "Predicted Cell Types (Dataset B)", 
       y = "Query Clusters (Dataset C)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
heatmap3

setwd("~/Google Drive/My Drive/Gensel lab/R studies/SingleR/")
save(singleR_result_BxC, file = "singleR_BxC_wilcox.RData")
save(df_long_BxC, file = "df_long_BxC_wilcox.RData")


##final plot --------------------------
common_limits <- range(c(df_long_AxB$value, df_long_AxC$value, df_long_BxC$value))
fill_scale <- scale_fill_gradient(low = "white", high = "#F4A3C0", limits = common_limits)
fill_scale <- scale_fill_gradient(low = "white", high = "#A1C9F4", limits = common_limits)

heatmap1 <- ggplot(df_long_AxB, aes(x = predicted_cell_types, y = query_clusters, fill = value)) +
  geom_tile(color = "white") +
  fill_scale +
  labs(title = "Dataset A\nvs Dataset B",
       x = "Dataset A", 
       y = "Dataset B",
       fill = "Number of cells\npredicted on\nreference clusters") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text (size = 15),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
heatmap1

heatmap2 <- ggplot(df_long_AxC, aes(x = predicted_cell_types, y = query_clusters, fill = value)) +
  geom_tile(color = "white") +
  fill_scale +
  labs(title = "Dataset A\nvs Dataset C",
       x = "Dataset A", 
       y = "Dataset C",
       fill = "Number of cells\npredicted on\nreference clusters") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text (size = 15),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
heatmap2

heatmap3 <- ggplot(df_long_BxC, aes(x = predicted_cell_types, y = query_clusters, fill = value)) +
  geom_tile(color = "white") +
  fill_scale +
  labs(title = "Dataset B\nvs Dataset C",
       x = "Dataset B", 
       y = "Dataset C",
       fill = "Number of cells\npredicted on\nreference clusters") +
  theme_minimal() +
  theme(plot.title = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text (size = 15),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
heatmap3

final_plot <- heatmap1 + heatmap2 + heatmap3 + 
  plot_layout(ncol = 3, guides = "collect", widths = c(1, 1, 1))
final_plot
setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 figures")
ggsave("SingleR_all_wilcox_2.png", dpi = 300, bg = "white")
