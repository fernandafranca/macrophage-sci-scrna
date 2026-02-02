# scripts/pipeline/04_singleR_cross_dataset.R
# Purpose: Cross-dataset macrophage cluster matching using SingleR (wilcox)
# Inputs (local-only cache from processing scripts):
#   - local_cache/salvador_A_macrophages_seurat.rds
#   - local_cache/milich_B_macrophage_monocyte_7dpi_seurat.rds
#   - local_cache/brennan_C_macrophage_monocyte_seurat.rds
# Outputs (repo):
#   - results/tables/singler/contingency_AxB.csv
#   - results/tables/singler/contingency_AxC.csv
#   - results/tables/singler/contingency_BxC.csv
#   - results/tables/singler/df_long_AxB.csv
#   - results/tables/singler/df_long_AxC.csv
#   - results/tables/singler/df_long_BxC.csv
#   - results/tables/singler/singler_labels_AxB.csv (per-cell labels, optional)
#   - results/tables/singler/singler_labels_AxC.csv
#   - results/tables/singler/singler_labels_BxC.csv
# Outputs (figures):
#   - figures/paper/singler_heatmaps.png

source("scripts/pipeline/00_config.R")

suppressPackageStartupMessages({
  library(SingleR)
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

# ---- Paths ----
out_dir <- file.path(DIR_TABLES, "singler")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fig_out <- file.path(DIR_FIGS, "paper")
dir.create(fig_out, showWarnings = FALSE, recursive = TRUE)

# ---- Load cached Seurat objects (local-only) ----
path_A <- file.path(DIR_LOCAL_CACHE, "salvador_A_macrophages_seurat.rds")                 # reference A
path_B <- file.path(DIR_LOCAL_CACHE, "milich_B_macro_seurat.rds")     # query B / reference B
path_C <- file.path(DIR_LOCAL_CACHE, "brennan_C_macro_seurat.rds")         # query C

stopifnot(file.exists(path_A), file.exists(path_B), file.exists(path_C))

mac_A <- readRDS(path_A)
mac_B <- readRDS(path_B)
mac_C <- readRDS(path_C)

# Ensure identities are set (SingleR labels = cluster identities)
Idents(mac_A) <- Idents(mac_A)
Idents(mac_B) <- Idents(mac_B)
Idents(mac_C) <- Idents(mac_C)

# -----------------------------------------------------------------------------
# Helper: run SingleR + make contingency + long df
# -----------------------------------------------------------------------------
run_singler_pair <- function(ref_obj, test_obj, pair_id,
                             ref_name = "ref", test_name = "test",
                             de_method = "wilcox") {
  
  ref_labels <- Idents(ref_obj)
  test_clusters <- as.factor(Idents(test_obj))
  
  # Run SingleR
  sr <- SingleR(
    test   = GetAssayData(test_obj, layer = "data"),
    ref    = GetAssayData(ref_obj,  layer = "data"),
    labels = ref_labels,
    de.method = de_method
  )
  
  # Per-cell predicted reference labels
  pred <- sr$labels
  
  # Contingency table: test cluster vs predicted reference label
  tab <- table(
    query_clusters = test_clusters,
    predicted_cell_types = pred
  )
  
  # Save contingency as CSV (wide)
  tab_df <- as.data.frame.matrix(tab) %>%
    rownames_to_column("query_cluster")
  
  write_csv(tab_df, file.path(out_dir, paste0("contingency_", pair_id, ".csv")))
  
  # Long format for plotting
  df_long <- as.data.frame(tab)
  
  # Robustly rename the count column to "value"
  # (it is usually "Freq", but not always)
  count_col <- setdiff(names(df_long), c("query_clusters", "predicted_cell_types"))[1]
  names(df_long)[names(df_long) == count_col] <- "value"
  
  df_long <- df_long %>%
    mutate(
      query_clusters = as.character(query_clusters),
      predicted_cell_types = as.character(predicted_cell_types),
      value = as.numeric(value),
      pair = pair_id,
      ref = ref_name,
      test = test_name
    )
  
  write_csv(df_long, file.path(out_dir, paste0("df_long_", pair_id, ".csv")))
  
  # Optional: save per-cell label assignments (can be big; still usually manageable)
  # If too large, comment these out.
  labels_df <- tibble(
    cell = colnames(test_obj),
    query_cluster = as.character(test_clusters),
    predicted_ref_cluster = pred
  )
  write_csv(labels_df, file.path(out_dir, paste0("singler_labels_", pair_id, ".csv")))
  
  list(singleR = sr, tab = tab, df_long = df_long)
}

# -----------------------------------------------------------------------------
# Run A vs B, A vs C, B vs C
# -----------------------------------------------------------------------------
res_AxB <- run_singler_pair(ref_obj = mac_A, test_obj = mac_B, pair_id = "AxB",
                            ref_name = "A", test_name = "B", de_method = "wilcox")

res_AxC <- run_singler_pair(ref_obj = mac_A, test_obj = mac_C, pair_id = "AxC",
                            ref_name = "A", test_name = "C", de_method = "wilcox")

res_BxC <- run_singler_pair(ref_obj = mac_B, test_obj = mac_C, pair_id = "BxC",
                            ref_name = "B", test_name = "C", de_method = "wilcox")

# -----------------------------------------------------------------------------
# Plot heatmaps with shared fill scale
# -----------------------------------------------------------------------------
df_long_AxB <- res_AxB$df_long
df_long_AxC <- res_AxC$df_long
df_long_BxC <- res_BxC$df_long

common_limits <- range(c(df_long_AxB$value, df_long_AxC$value, df_long_BxC$value))
fill_scale <- scale_fill_gradient(low = "white", high = "#F4A3C0", limits = common_limits)

make_heatmap <- function(df_long, title, xlab, ylab) {
  ggplot(df_long, aes(x = predicted_cell_types, y = query_clusters, fill = value)) +
    geom_tile(color = "white") +
    fill_scale +
    labs(
      title = title,
      x = xlab,
      y = ylab,
      fill = "Number of cells\npredicted on\nreference clusters"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18),
      axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 13),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

p1 <- make_heatmap(df_long_AxB, "Dataset A\nvs Dataset B", xlab = "Dataset A", ylab = "Dataset B")
p2 <- make_heatmap(df_long_AxC, "Dataset A\nvs Dataset C", xlab = "Dataset A", ylab = "Dataset C")
p3 <- make_heatmap(df_long_BxC, "Dataset B\nvs Dataset C", xlab = "Dataset B", ylab = "Dataset C")

final_plot <- p1 + p2 + p3 + plot_layout(ncol = 3, guides = "collect")

ggsave(
  filename = file.path(fig_out, "singler_heatmaps.png"),
  plot = final_plot,
  width = 18,
  height = 6,
  dpi = 300,
  bg = "white"
)

message("✅ SingleR cross-dataset matching complete.")
