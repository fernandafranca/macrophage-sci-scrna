library(tidyverse)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(SeuratDisk)
library(HDF5Array)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")

#get data
mtx_obj <- ReadMtx(mtx = "data/GSM5904827_B3_matrix.mtx.gz",
                   features = "data/GSM5904827_B3_features.tsv.gz",
                   cells = "data/GSM5904827_B3_barcodes.tsv.gz")
seurat_mtx <- CreateSeuratObject(counts = mtx_obj, min.cells = 3, min.features = 200)
#15313  7439

#QC and normalization
seurat_mtx[["percent.mt"]] <- PercentageFeatureSet(seurat_mtx, pattern = "^mt-")
seurat_mtx <- subset(seurat_mtx, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & 
                       percent.mt < 20 & nCount_RNA > 1500 & nCount_RNA < 50000)
# 15313  6478
# 15313 features across 6458 sample for nfeatures >500

seurat_mtx <- NormalizeData(seurat_mtx)
seurat_mtx <- FindVariableFeatures(seurat_mtx, selection.method = "vst")
all.genes <- rownames(seurat_mtx)
seurat_mtx <- ScaleData(seurat_mtx, features = all.genes)
seurat_mtx <- RunPCA(seurat_mtx, features = VariableFeatures(object = seurat_mtx))
ElbowPlot(seurat_mtx)

# Clustering 
seurat_mtx <- FindNeighbors(seurat_mtx, dims = 1:15)
seurat_mtx <- FindClusters(seurat_mtx, resolution = 0.1) #confirm with someone about this resolution
DimPlot(seurat_mtx, group.by = "RNA_snn_res.0.1", label = TRUE)
Idents(seurat_mtx) <- "RNA_snn_res.0.1"

#tSNE
seurat_mtx <- RunTSNE(seurat_mtx, dims = 1:15)
pastel_colors <- c("#A1C9F4", "#FFB482", "#8DE5A1", "#F4A3C0", "#D0BBFF", 
                   "#FDE051", "#FB8072", "#9DDECB", "#C6A78B", "#BAB0AC")

TSNEPlot(seurat_mtx, raster = F, cols =pastel_colors)
 
#umap
seurat_mtx <- RunUMAP(seurat_mtx, dims = 1:15)
DimPlot(seurat_mtx, reduction = "umap", label = TRUE)

## get cluster markers
cluster1.markers <- FindAllMarkers(seurat_mtx, logfc.threshold = .25, min.pct= 0.3, only.pos = T)
cluster1.markers <- cluster1.markers[which(cluster1.markers$p_val_adj < 0.05), ]
#use canonical markers to identify clusters
#Brennan markers
DotPlot(seurat_mtx, features = c("P2ry12","Siglech","Tmem119","Gpnmb","Spp1","Fabp5","Ly6c1","Cldn5","Itm2a","H2-Aa","H2-Eb1","H2-Ab1","Slc1a2","Atp1a2","Atp1b2","Ccl5","Ms4a4b","Nkg7","Dbi","Nnat","Mt3","Ly6d","Cd79a","Cd79b","Stmn1","Hmgb2","Top2a","S100a9","S100a8","Retnlg","Hba-a2","Hbb-bs","Alas2","Olig1","Cspg5","Vcan","Apod","Dcn","Vtn","Mgp","Igfbp6","Gsn"),
        cols =  c("#704D9E", "#EB7999")) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 17, angle = 45, hjust = 1), axis.text.y = element_text (size = 15), text = element_text(size = 10))
ggsave("brennan_dotplot_named_colors_size.png", width = 6, height = 10)

#Salvador markers
DotPlot(seurat_mtx, features = c("Csf1r","Cd14","Cd68","Ptprc","Tmem119","Sall1","P2ry12","Ms4a7","Adgre1","Ms4a6c","Ly6c2","Ly6g","S100a9","H2-Ab1","Cd74","Itgax","Cd3e","Ms4a4b","Klrb1c","Ly6d","Igkc","Cd19","Cldn5","Col4a1","Pdgfrb","Rgs5")) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 12))
ggsave("dotplot_named_salvador_markers.png", width = 6, height = 10)


 #add proposed annotation
new.clusters <- c('Microglia','Macrophages', 'Monocytes', 'Endothelial Cells', 'Intermediate Progenitors','Ependymal cells','Neutrophils', 'Oligodendrocytes', "Miscelaneous", "Astrocytes")
names(new.clusters) <- levels(seurat_mtx)
seurat_mtx <- RenameIdents(seurat_mtx, new.clusters)

#save new identity and visualize
seurat_mtx$celltype <- seurat_mtx@active.ident
TSNEPlot(seurat_mtx, raster = F)

save(seurat_mtx, file = "seurat_mtx_gse196928.RData")

## ------- re-cluster macs -----------------------
macro <- subset(seurat_mtx, ident = c("Macrophages", "Monocytes"))
macro <- ScaleData(macro, features = row.names(macro),
                  vars.to.regress = c("nFeature_RNA","nCount_RNA","percent_mito"))
macro <- FindVariableFeatures(macro, selection.method = "vst")
macro <- RunPCA(macro, npcs = 50, features = VariableFeatures(macro))
ElbowPlot(macro, 50)

macro <- FindNeighbors(macro, dims = 1:30, reduction = "pca")
macro <- FindClusters(macro, resolution = 0.4)
DimPlot(macro, group.by = "RNA_snn_res.0.4", label = TRUE)
macro <- RunTSNE(macro, dims = 1:30, verbose = T)

TSNEPlot(macro, cols =pastel_colors)

## characterize clusters
clustmarks <- FindAllMarkers(macro, logfc.threshold = .25, min.pct= 0.3, only.pos = T)
clustmarks <- clustmarks[which(clustmarks$p_val_adj < 0.05), ]

#top 10 by adj p value
clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10
DoHeatmap(macro, features = top10$gene, group.colors = pastel_colors)+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
write.csv(top10, file = "top10_brennan_adjpval_2.csv")

#top 10 by logfc
clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10_2
DoHeatmap(macro, features = top10_2$gene, group.colors = paletteer_c("grDevices::Sunset", 8)) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) 
write.csv(top10_2, file = "top10_brennan_logfc_2.csv")
save(macro, file = "macro_brennan.RData")

###reviewer's suggestions ---------
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
load("macro_brennan.RData")

#Markers plots
VlnPlot(macro, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl","Sparc","Cd81","Gpr34","Serpine2","Gpnmb","Cdo1","Hpse","Ccr2","Plac8"), pt.size = 0) 
DotPlot(macro, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl","Sparc","Cd81","Gpr34","Serpine2","Gpnmb","Cdo1","Hpse","Ccr2","Plac8")) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text (size = 20), text = element_text(size = 20))
ggsave("brennan_dotplot_markers.png", width = 6, height = 10)
FeaturePlot(macro, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl","Sparc","Cd81","Gpr34","Serpine2","Gpnmb","Cdo1","Hpse","Ccr2","Plac8"), reduction = "tsne")

#number and proportions of cells per cluster
cluster_counts <- table(macro@active.ident)
cluster_proportions <- prop.table(cluster_counts)
cluster_summary_brennan <- data.frame(
  Cluster = names(cluster_counts),
  Cell_Count = as.numeric(cluster_counts),
  Proportion = as.numeric(cluster_proportions)
)
save(cluster_summary_brennan, file = "cluster_summary_brennan.RData")
write.csv(cluster_summary_brennan, file = "cluster_summary_brennan.csv")

ggplot(cluster_summary_brennan, aes(x = Cluster, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Cells in Each Cluster", y = "Proportion", x = "Cluster") +
  theme_minimal()


