# Dataset A — Salvador et al. (GSE205037), Young 7dpi
# Purpose: QC → clustering → coarse cell-type labels → macrophage subset reclustering
# Outputs:
#   - results/tables/cluster_summary_salvador.csv
#   - results/tables/top10_salvador_macrophage_clusters.csv
#   - figures/salvador_A/* (representative plots)
#
# Notes:
# - Raw data and large objects are stored locally (local_data/, local_cache/) and not tracked by git.

source("scripts/pipeline/00_config.R")

REGRESS_MITO <- FALSE   # paper-matching
set.seed(1234)


# ---- Libraries used in this script  ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
})

dataset_id <- "salvador_A"

# ---- Local inputs ----
mtx_dir <- file.path(DIR_LOCAL_DATA, "GSE205037_RAW_young", "Young.7d")

mtx_file      <- file.path(mtx_dir, "matrix.mtx.gz")
features_file <- file.path(mtx_dir, "features.tsv.gz")
barcodes_file <- file.path(mtx_dir, "barcodes.tsv.gz")

stopifnot(file.exists(mtx_file), file.exists(features_file), file.exists(barcodes_file))

# ---- Output folders ----
fig_dir <- file.path(DIR_FIGS, "salvador_A")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1) Load raw matrix + create Seurat object
# =============================================================================
mtx_obj <- ReadMtx(mtx = mtx_file, features = features_file, cells = barcodes_file)

seu <- CreateSeuratObject(counts = mtx_obj, min.cells = 3, min.features = 200)
message("Loaded Salvador 7dpi: ", nrow(seu), " genes x ", ncol(seu), " cells")

# =============================================================================
# 2) QC filtering 
# =============================================================================
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

seu <- subset(
  seu,
  subset = nFeature_RNA > 500 &
    nFeature_RNA < 7500 &
    percent.mt < 20 &
    nCount_RNA > 1500 &
    nCount_RNA < 50000
)

message("After QC: ", nrow(seu), " genes x ", ncol(seu), " cells")

# =============================================================================
# 3) Normalize + HVGs + scale + PCA
# =============================================================================
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst")

seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(seu))

DimHeatmap(seu, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(seu)


# =============================================================================
# 4) Clustering + embeddings
# =============================================================================
dims_use <- 1:13

seu <- FindNeighbors(seu, dims = dims_use)
seu <- FindClusters(seu, resolution = 0.05)
Idents(seu) <- "RNA_snn_res.0.05"

seu <- RunTSNE(seu, dims = dims_use)
seu <- RunUMAP(seu, dims = dims_use)

# =============================================================================
# 5) Cluster markers (all cells) + coarse annotation
# =============================================================================
all_markers <- FindAllMarkers(seu, logfc.threshold = 0.25, min.pct = 0.3, only.pos = TRUE) %>%
  filter(p_val_adj < 0.05)

# write per-cluster markers (top by avg_log2FC)
dir.create(file.path(DIR_TABLES, "salvador_A"), showWarnings = FALSE, recursive = TRUE)

for (cl in sort(unique(all_markers$cluster))) {
  tmp <- all_markers %>%
    filter(cluster == cl) %>%
    arrange(desc(avg_log2FC))
  write.csv(tmp, file.path(DIR_TABLES, "salvador_A", paste0("markers_allcells_cluster_", cl, ".csv")),
            row.names = FALSE)
}

# Canonical marker dotplot for broad annotation
p_dot <- DotPlot(
  seu,
  features = c(
    "Csf1r","Cd14","Cd68","Ptprc",
    "Tmem119","Sall1","P2ry12",
    "Ms4a7","Adgre1","Ms4a6c",
    "Ly6c2","Ly6g","S100a9",
    "H2-Ab1","Cd74","Itgax",
    "Cd3e","Ms4a4b","Klrb1c",
    "Ly6d","Igkc","Cd19",
    "Cldn5","Col4a1","Pdgfrb","Rgs5"
  )
) + coord_flip()


# Coarse labels
new_clusters <- c("Microglia", "Macrophages", "Neutrophils", "Lymphocytes", "Stromal/ECs")
names(new_clusters) <- levels(seu)
seu <- RenameIdents(seu, new_clusters)
seu$celltype <- Idents(seu)

p_tsne_celltype <- TSNEPlot(seu, raster = FALSE) +
  ggtitle("Salvador (A) — coarse cell types (tSNE)")

# =============================================================================
# 6) Subset macrophages and recluster
# =============================================================================
myel <- subset(seu, ident = "Macrophages")

vars_regress <- c("nFeature_RNA", "nCount_RNA")
if (REGRESS_MITO) vars_regress <- c(vars_regress, "percent.mt")


myel <- ScaleData(
  myel,
  features = rownames(myel),
  vars.to.regress = vars_regress
)

myel <- FindVariableFeatures(myel, selection.method = "vst")
myel <- RunPCA(myel, npcs = 50, features = VariableFeatures(myel))

ElbowPlot(myel, ndims = 50)

myel_dims <- 1:40
myel <- FindNeighbors(myel, dims = myel_dims, reduction = "pca")
myel <- FindClusters(myel, resolution = 0.4)
myel <- RunTSNE(myel, dims = myel_dims)

p_tsne_myel <- TSNEPlot(myel) +
  ggtitle("Salvador (A) — macrophages reclustered (tSNE)")


# =============================================================================
# 7) Macrophage cluster markers + top10 signatures 
# =============================================================================
mac_markers <- FindAllMarkers(myel, logfc.threshold = 0.25, min.pct = 0.3, only.pos = TRUE) %>%
  filter(p_val_adj < 0.05)

# Top10: adjP<0.05 AND avg_log2FC>1, then first 10 per cluster
top10 <- mac_markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(top10, file.path(DIR_TABLES, "top10_salvador_macrophage_clusters.csv"), row.names = FALSE)

p_hm1 <- DoHeatmap(myel, features = unique(top10$gene)) +
  ggtitle("Salvador (A) macrophages — Top markers (avg_log2FC>1)")

# Top10 by logFC: sort by avg_log2FC within cluster
top10_logfc <- mac_markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 10) %>%
  ungroup()

write.csv(top10_logfc, file.path(DIR_TABLES, "top10_logfc_salvador_macrophage_clusters.csv"), row.names = FALSE)

# =============================================================================
# 8) Cluster counts + proportions (macrophage reclustering)
# =============================================================================
cluster_counts <- table(Idents(myel))
cluster_props  <- prop.table(cluster_counts)

cluster_summary_salvador <- tibble(
  Cluster = names(cluster_counts),
  Cell_Count = as.numeric(cluster_counts),
  Proportion = as.numeric(cluster_props)
)

write_csv(cluster_summary_salvador, file.path(DIR_TABLES, "cluster_summary_salvador.csv"))

p_bar <- ggplot(cluster_summary_salvador, aes(x = Cluster, y = Proportion, fill = Cluster)) +
  geom_col() +
  labs(title = "Salvador (A) macrophages: cluster proportions", y = "Proportion", x = "Cluster") +
  theme_minimal()


# =============================================================================
# 9) Save big objects locally 
# =============================================================================
saveRDS(seu,  file.path(DIR_LOCAL_CACHE, "salvador_A_allcells_seurat.rds"))
message("REGRESS_MITO = ", REGRESS_MITO)
message("nCells = ", ncol(myel))
message("nClusters = ", length(levels(Idents(myel))))
print(sort(table(Idents(myel)), decreasing = TRUE))

saveRDS(myel, file = file.path(DIR_LOCAL_CACHE,"salvador_A_macrophages_seurat.rds"))



