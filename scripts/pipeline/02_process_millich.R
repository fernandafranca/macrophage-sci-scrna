# Dataset B — Milich et al. (GSE162610), 7dpi (myeloid subset from authors)
# Purpose: subset 7dpi → macrophage/monocyte subset → recluster → markers → summary tables
# Inputs (local-only): local_data/GSE162610_Milich/myeloid-003.rds
# Outputs (repo):
#   - results/tables/cluster_summary_milich.csv
#   - results/tables/top10_milich_macrophage_clusters.csv
#   - results/tables/top10_logfc_milich_macrophage_clusters.csv
#   - results/tables/milich_B/markers_cluster_<id>.csv  (optional)
# Outputs (local-only cache, not tracked):
#   - local_cache/milich_B_myeloid_7dpi_seurat.rds
#   - local_cache/milich_B_macrophage_monocyte_7dpi_seurat.rds

source("scripts/pipeline/00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
})

dataset_id <- "milich_B"

# -----------------------------------------------------------------------------
# 1) Load myeloid object and subset 7dpi
# -----------------------------------------------------------------------------
input_rds <- file.path(DIR_LOCAL_DATA, "GSE162610_Milich", "myeloid-003.rds")
stopifnot(file.exists(input_rds))

myel <- readRDS(input_rds)

# Subset to 7dpi samples
Idents(myel) <- "sample_id"
myel <- subset(myel, idents = c("7dpi_sample1", "7dpi_sample2"))

# Set identity to celltype for downstream subsetting
Idents(myel) <- "celltype"
ElbowPlot(myel, 20)

# Compute tSNE
myel <- RunTSNE(myel, dims = 1:15, verbose = T)
TSNEPlot(myel)


# Cache object locally (NOT tracked by git)
saveRDS(myel, file.path(DIR_LOCAL_CACHE, "milich_B_myeloid_seurat.rds"))

# -----------------------------------------------------------------------------
# 2) Subset macrophages/monocytes and recluster
# -----------------------------------------------------------------------------
# These label names must match your object’s "celltype" levels
macro <- subset(myel, idents = c("Macrophage", "Monocyte"))

macro <- NormalizeData(macro)
macro <- FindVariableFeatures(macro, selection.method = "vst")

# Regress technical covariates 
macro <- ScaleData(
  macro,
  features = rownames(macro),
  vars.to.regress = c("nFeature_RNA", "nCount_RNA")
)

macro <- RunPCA(macro, npcs = 50, features = VariableFeatures(macro))
ElbowPlot(macro, ndims = 50)

macro <- FindNeighbors(macro, dims = 1:30, reduction = "pca")
macro <- FindClusters(macro, resolution = 0.31)

macro <- RunTSNE(macro, dims = 1:40, verbose = TRUE)
TSNEPlot(macro)

# Cache macrophage/monocyte object locally
saveRDS(macro, file.path(DIR_LOCAL_CACHE, "milich_B_macro_seurat.rds"))

# -----------------------------------------------------------------------------
# 3) Markers + top10 signatures
# -----------------------------------------------------------------------------
clustmarks <- FindAllMarkers(macro, logfc.threshold = 0.25, min.pct = 0.3, only.pos = TRUE) %>%
  filter(p_val_adj < 0.05)

# Write full marker tables per cluster 
dir.create(file.path(DIR_TABLES, "milich_B"), showWarnings = FALSE, recursive = TRUE)

for (cl in sort(unique(clustmarks$cluster))) {
  tmp <- clustmarks %>%
    filter(cluster == cl) %>%
    arrange(desc(avg_log2FC))
  write.csv(tmp,
            file.path(DIR_TABLES, "milich_B", paste0("markers_cluster_", cl, ".csv")),
            row.names = FALSE)
}

# Top10 by adj p (avg_log2FC > 1, then first 10 per cluster)
top10 <- clustmarks %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(top10, file.path(DIR_TABLES, "top10_milich_macrophage_clusters.csv"), row.names = FALSE)

# Top10 by logFC (sorted)
top10_logfc <- clustmarks %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(top10_logfc, file.path(DIR_TABLES, "top10_logfc_milich_macrophage_clusters.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 4) Cluster counts + proportions
# -----------------------------------------------------------------------------
cluster_counts <- table(Idents(macro))
cluster_props  <- prop.table(cluster_counts)

cluster_summary_milich <- tibble(
  Cluster = names(cluster_counts),
  Cell_Count = as.numeric(cluster_counts),
  Proportion = as.numeric(cluster_props)
)

write_csv(cluster_summary_milich, file.path(DIR_TABLES, "cluster_summary_milich.csv"))


