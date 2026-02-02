# scripts/figures/01_paper_figures.R
# Purpose: Final, styled figures for the macrophage scRNA paper 
# Inputs: locally cached Seurat objects 
# Outputs: figures/paper/*.png

source("scripts/pipeline/00_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(Seurat)
  library(ggh4x)
})

# -------------------- Styling --------------------
pastel_colors <- c("#A1C9F4", "#F4A3C0", "#FFB482", "#8DE5A1", "#D0BBFF",
                   "#FDE051", "#FB8072", "#9DDECB", "#C6A78B", "#BAB0AC")
gradient_colors <- c("#A1C9F4", "white", "#F4A3C0")

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

# Where figures will be saved 
fig_out <- file.path(DIR_FIGS, "paper")
dir.create(fig_out, showWarnings = FALSE, recursive = TRUE)

# -------------------- Salvador (A) --------------------
# (These were created by scripts/pipeline/01_process_salvador.R)
salvador_allcells_path <- file.path(DIR_LOCAL_CACHE, "salvador_A_allcells_seurat.rds")
stopifnot(file.exists(salvador_allcells_path))

seurat_mtx <- readRDS(salvador_allcells_path)

# All cells dotplot
features_salvador <- c(
  "Csf1r","Cd14","Cd68","Ptprc",
  "Tmem119","Sall1","P2ry12",
  "Ms4a7","Adgre1","Ms4a6c",
  "Ly6c2","Ly6g","S100a9",
  "H2-Ab1","Cd74","Itgax",
  "Cd3e","Ms4a4b","Klrb1c",
  "Ly6d","Igkc","Cd19",
  "Cldn5","Col4a1","Pdgfrb","Rgs5"
)

p_dot_salvador <- DotPlot(seurat_mtx, features = features_salvador) +
  scale_color_gradient2(low = gradient_colors[1], mid = gradient_colors[2], high = gradient_colors[3]) +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    text = element_text(size = 15),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )

ggsave(
  filename = file.path(fig_out, "dotplot_salvador_allcells.png"),
  plot = p_dot_salvador,
  width = 8,
  height = 10,
  bg = "white",
  dpi = 300
)

#All cells TSNE
p_tsne_salvador <- TSNEPlot(seurat_mtx, raster = F, cols =pastel_colors)+
  guides (x = axis, y = axis, color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3)))+
  labs(x = "tSNE 1", y = "tSNE 2", title = "Spinal cord myeloid cells")+
  theme(axis.title = element_text(hjust = 0),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 20, face = "plain", hjust = 0.5))
ggsave(
  filename = file.path(fig_out, "tsne_salvador_all.png"),
  plot = p_tsne_salvador,
  width = 8,
  height = 10,
  bg = "white",
  dpi = 300
)

#Macrophages tSNE
salvador_myel_path <- file.path(DIR_LOCAL_CACHE, "salvador_A_macrophages_seurat.rds")
myel <- readRDS(salvador_myel_path)

p_tsne_myel_salvador <- TSNEPlot(myel, cols = pastel_colors)+
  guides (x = axis, y = axis, color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 3)))+
  labs(x = "tSNE 1", y = "tSNE 2", title = "Macrophages")+
  theme(axis.title = element_text(hjust = 0),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(size = 20, face = "plain", hjust = 0.5))

ggsave(
  filename = file.path(fig_out, "tsne_salvador_macrophages.png"),
  plot = p_tsne_myel_salvador,
  width = 8,
  height = 10,
  bg = "white",
  dpi = 300
)

##Heatmaps
clustmarks <- FindAllMarkers(myel, logfc.threshold = .25, min.pct= 0.3, only.pos = T)
clustmarks <- clustmarks[which(clustmarks$p_val_adj < 0.05), ]

genes_to_highlight <- unique(c("Fabp4","Lgals1","Lgals3","Spp1","Cd81","Gpr34","Serpine2","Sparc","Ahnak2","Cdo1","Gpnmb","Hpse","Ccr2","Ifitm6","Plac8","Mrc1","Cd163","Cd36","Lpl","Mfge8","Spp1","Gpr84","Sparc","Cdo1","Plpp3","Ifitm6","Mcemp1","Plac8","Apoc1","Cd163","Clec4a1","Mrc1"))

clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10

heatmap_salvador <- DoHeatmap(myel, features = top10$gene,  group.colors = pastel_colors) +
  scale_fill_gradientn(colors =  rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) 
heatmap_salvador + theme(axis.text.y = element_text(color = ifelse(levels(p1$data$Feature) %in% genes_to_highlight, "red", "black")))

ggsave(
  filename = file.path(fig_out, "heatmap_salvador_pvalue.png"),
  width = 8,
  height = 10,
  bg = "white",
  dpi = 300
)

#top10 by logfc
clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10_2

heatmap_salvador_logfc <- DoHeatmap(myel, features = top10_2$gene, group.colors = pastel_colors)+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
heatmap_salvador_logfc + theme(axis.text.y = element_text(color = ifelse(levels(p2$data$Feature) %in% genes_to_highlight, "red", "black")))

ggsave(
  filename = file.path(fig_out, "heatmap_salvador_logfc.png"),
  width = 8,
  height = 10,
  bg = "white",
  dpi = 300
)

# -------------------- Milich (B) --------------------
# (These were created by scripts/pipeline/01_process_milich.R)
milich_allcells_path <- file.path(DIR_LOCAL_CACHE, "milich_B_myeloid_seurat.rds")
milich_macro_path    <- file.path(DIR_LOCAL_CACHE, "milich_B_macro_seurat.rds")

stopifnot(file.exists(milich_allcells_path))
stopifnot(file.exists(milich_macro_path))

myel_milich <- readRDS(milich_allcells_path)
macro_milich <- readRDS(milich_macro_path)

# All cells dotplot
features_milich <- c(
  "Gpr84","Ptgs1","P2ry12","Gpr34","Lag3","Siglech","Cst7",
  "Ly6c2","Ccr2","F10","Plac8","Thbs1","Ms4a7",
  "Pf4","Fabp4","Gpnmb",
  "Top2a","Mki67","Cdk1","Ccnb2",
  "Cd74","S100a9","Mmp9","Ly6g","Cd177","Ltf"
)

p_dot_milich <- DotPlot(myel_milich, features = features_milich) +
  scale_color_gradient2(low = gradient_colors[1], mid = gradient_colors[2], high = gradient_colors[3]) +
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    text = element_text(size = 15)
  )

ggsave(file.path(fig_out, "dotplot_milich_allcells.png"),
       p_dot_milich, width = 8, height = 10, dpi = 300)

# All cells tSNE
p_tsne_milich <- TSNEPlot(myel_milich, cols = pastel_colors) +
  guides(x = axis, y = axis) +
  labs(title = "Spinal cord myeloid cells") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

ggsave(file.path(fig_out, "tsne_milich_all.png"),
       p_tsne_milich, width = 8, height = 6, dpi = 300)

# Macrophages tSNE
p_tsne_milich_macro <- TSNEPlot(macro_milich, cols = pastel_colors) +
  guides(x = axis, y = axis) +
  labs(title = "Macrophages") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

ggsave(file.path(fig_out, "tsne_milich_macrophages.png"),
       p_tsne_milich_macro, width = 8, height = 6, dpi = 300)

# Heatmaps
clustmarks_milich <- FindAllMarkers(macro_milich, only.pos = TRUE) %>%
  filter(p_val_adj < 0.05)

top10_milich <- clustmarks_milich %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

heatmap_milich <- DoHeatmap(macro_milich, features = top10_milich$gene,
                            group.colors = pastel_colors)

ggsave(file.path(fig_out, "heatmap_milich_logfc.png"),
       heatmap_milich, width = 8, height = 10, dpi = 300)
# -------------------- Brennan (C) --------------------
brennan_allcells_path <- file.path(DIR_LOCAL_CACHE, "brennan_C_allcells_seurat.rds")
brennan_macro_path    <- file.path(DIR_LOCAL_CACHE, "brennan_C_macro_seurat.rds")

stopifnot(file.exists(brennan_allcells_path))
stopifnot(file.exists(brennan_macro_path))

seurat_brennan <- readRDS(brennan_allcells_path)
macro_brennan  <- readRDS(brennan_macro_path)

# All cells dotplot
features_brennan <- c(
  "P2ry12","Siglech","Tmem119","Gpnmb","Spp1","Fabp5",
  "Ly6c1","Cldn5","Itm2a",
  "H2-Aa","H2-Eb1","H2-Ab1",
  "Slc1a2","Atp1a2","Atp1b2",
  "Ccl5","Ms4a4b","Nkg7",
  "Stmn1","Hmgb2","Top2a",
  "S100a9","S100a8","Retnlg",
  "Vcan","Apod","Dcn","Vtn","Mgp"
)

p_dot_brennan <- DotPlot(seurat_brennan, features = features_brennan) +
  scale_color_gradient2(low = gradient_colors[1], mid = gradient_colors[2], high = gradient_colors[3]) +
  coord_flip()

ggsave(file.path(fig_out, "dotplot_brennan_allcells.png"),
       p_dot_brennan, width = 8, height = 10, dpi = 300)

# All cells tSNE
p_tsne_brennan <- TSNEPlot(seurat_brennan, cols = pastel_colors) +
  guides(x = axis, y = axis) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

ggsave(file.path(fig_out, "tsne_brennan_all.png"),
       p_tsne_brennan, width = 8, height = 6, dpi = 300)

# Macrophages tSNE
p_tsne_brennan_macro <- TSNEPlot(macro_brennan, cols = pastel_colors) +
  guides(x = axis, y = axis) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom")

ggsave(file.path(fig_out, "tsne_brennan_macrophages.png"),
       p_tsne_brennan_macro, width = 8, height = 6, dpi = 300)

# Heatmaps
clustmarks_brennan <- FindAllMarkers(macro_brennan, only.pos = TRUE) %>%
  filter(p_val_adj < 0.05)

top10_brennan <- clustmarks_brennan %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

heatmap_brennan <- DoHeatmap(macro_brennan, features = top10_brennan$gene,
                             group.colors = pastel_colors)

ggsave(file.path(fig_out, "heatmap_brennan_logfc.png"),
       heatmap_brennan, width = 8, height = 10, dpi = 300)

