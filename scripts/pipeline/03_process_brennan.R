# Dataset C — Brennan et al. (GSE196928), sample B3 (7dpi)
# Purpose: QC → clustering → coarse annotation → macrophage/monocyte subset reclustering → markers + tables
# Inputs (local-only):
#   local_data/GSE196928/data/GSM5904827_B3_matrix.mtx.gz
#   local_data/GSE196928/data/GSM5904827_B3_features.tsv.gz
#   local_data/GSE196928/data/GSM5904827_B3_barcodes.tsv.gz
# Outputs (repo):
#   - results/tables/cluster_summary_brennan.csv
#   - results/tables/top10_brennan_macrophage_clusters.csv
#   - results/tables/top10_logfc_brennan_macrophage_clusters.csv
#   - results/tables/brennan_C/markers_allcells_cluster_<id>.csv (optional)
# Outputs (local-only cache):
#   - local_cache/brennan_C_allcells_seurat.rds
#   - local_cache/brennan_C_macrophage_monocyte_seurat.rds

source("scripts/pipeline/00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
})

dataset_id <- "brennan_C"

# -----------------------------------------------------------------------------
# 1) Load matrix and create Seurat object
# -----------------------------------------------------------------------------
mtx_dir <- file.path(DIR_LOCAL_DATA, "GSE196928", "data")

mtx_file      <- file.path(mtx_dir, "GSM5904827_B3_matrix.mtx.gz")
features_file <- file.path(mtx_dir, "GSM5904827_B3_features.tsv.gz")
barcodes_file <- file.path(mtx_dir, "GSM5904827_B3_barcodes.tsv.gz")

stopifnot(file.exists(mtx_file), file.exists(features_file), file.exists(barcodes_file))

mtx_obj <- ReadMtx(mtx = mtx_file, features = features_file, cells = barcodes_file)

seu <- CreateSeuratObject(counts = mtx_obj, min.cells = 3, min.features = 200)
message("Loaded Brennan B3: ", nrow(seu), " genes x ", ncol(seu), " cells")

# -----------------------------------------------------------------------------
# 2) QC filtering 
# -----------------------------------------------------------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

seu <- subset(
  seu,
  subset =
    nFeature_RNA > 500 & nFeature_RNA < 7500 &
    percent.mt < 20 &
    nCount_RNA > 1500 & nCount_RNA < 50000
)

message("After QC: ", nrow(seu), " genes x ", ncol(seu), " cells")

# -----------------------------------------------------------------------------
# 3) Normalize + HVGs + scale + PCA
# -----------------------------------------------------------------------------
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst")

seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(seu))
ElbowPlot(seu)

# -----------------------------------------------------------------------------
# 4) Clustering + embeddings 
# -----------------------------------------------------------------------------
dims_use <- 1:15

seu <- FindNeighbors(seu, dims = dims_use)
seu <- FindClusters(seu, resolution = 0.1)
Idents(seu) <- "RNA_snn_res.0.1"

# Compute embeddings for later figure script
set.seed(123)
seu <- RunTSNE(seu, dims = dims_use)
TSNEPlot(seu)

# -----------------------------------------------------------------------------
# 5) Cluster markers (all cells) — save per cluster 
# -----------------------------------------------------------------------------
all_markers <- FindAllMarkers(seu, logfc.threshold = 0.25, min.pct = 0.3, only.pos = TRUE) %>%
  filter(p_val_adj < 0.05)

dir.create(file.path(DIR_TABLES, "brennan_C"), showWarnings = FALSE, recursive = TRUE)

for (cl in sort(unique(all_markers$cluster))) {
  tmp <- all_markers %>%
    filter(cluster == cl) %>%
    arrange(desc(avg_log2FC))
  write.csv(tmp,
            file.path(DIR_TABLES, "brennan_C", paste0("markers_allcells_cluster_", cl, ".csv")),
            row.names = FALSE)
}

# -----------------------------------------------------------------------------
# 6) Coarse cell-type annotation 
# -----------------------------------------------------------------------------
new_clusters <- c(
  "Microglia", "Macrophages", "Monocytes", "Endothelial Cells",
  "Intermediate Progenitors", "Ependymal cells", "Neutrophils",
  "Oligodendrocytes", "Miscelaneous", "Astrocytes"
)

names(new_clusters) <- levels(seu)
seu <- RenameIdents(seu, new_clusters)
seu$celltype <- Idents(seu)

# Cache all-cells object locally
saveRDS(seu, file.path(DIR_LOCAL_CACHE, "brennan_C_allcells_seurat.rds"))

# -----------------------------------------------------------------------------
# 7) Subset macrophages + monocytes and recluster
# -----------------------------------------------------------------------------
macro <- subset(seu, ident = c("Macrophages", "Monocytes"))
macro <- NormalizeData(macro)
macro <- FindVariableFeatures(macro, selection.method = "vst")

macro <- ScaleData(
  macro,
  features = rownames(macro),
  vars.to.regress = c("nFeature_RNA", "nCount_RNA")
)

macro <- RunPCA(macro, npcs = 50, features = VariableFeatures(macro))
ElbowPlot(macro, ndims = 50)

macro <- FindNeighbors(macro, dims = 1:30, reduction = "pca")
macro <- FindClusters(macro, resolution = 0.4)

# Compute embeddings for later figure script
set.seed(123)
macro <- RunTSNE(macro, dims = 1:30, verbose = TRUE)
TSNEPlot(macro)

# Cache macrophage+monocyte object locally
saveRDS(macro, file.path(DIR_LOCAL_CACHE, "brennan_C_macro_seurat.rds"))

# -----------------------------------------------------------------------------
# 8) Macrophage cluster markers + top10 signatures (tables only)
# -----------------------------------------------------------------------------
clustmarks <- FindAllMarkers(macro, logfc.threshold = 0.25, min.pct = 0.3, only.pos = TRUE) %>%
  filter(p_val_adj < 0.05)

# Top10 by adj p (with avg_log2FC > 1) 
top10 <- clustmarks %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(top10, file.path(DIR_TABLES, "top10_brennan_macrophage_clusters.csv"), row.names = FALSE)

# Top10 by logFC — as legacy
top10_logfc <- clustmarks %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(top10_logfc, file.path(DIR_TABLES, "top10_logfc_brennan_macrophage_clusters.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 9) Cluster counts + proportions (macrophage+monocyte reclustering)
# -----------------------------------------------------------------------------
cluster_counts <- table(Idents(macro))
cluster_props  <- prop.table(cluster_counts)

cluster_summary_brennan <- tibble(
  Cluster = names(cluster_counts),
  Cell_Count = as.numeric(cluster_counts),
  Proportion = as.numeric(cluster_props)
)

write_csv(cluster_summary_brennan, file.path(DIR_TABLES, "cluster_summary_brennan.csv"))

message("✅ Brennan dataset processing complete.")
