library(tidyverse)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(SeuratDisk)
library(HDF5Array)
library(paletteer)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/")
mtx_obj <- ReadMtx(mtx = "GSE205037_RAW_young/Young.7d/matrix.mtx.gz",
                   features = "GSE205037_RAW_young/Young.7d/features.tsv.gz",
                   cells = "GSE205037_RAW_young/Young.7d/barcodes.tsv.gz")
seurat_mtx <- CreateSeuratObject(counts = mtx_obj, min.cells = 3, min.features = 200)
#17858 features across 28228 sample

# 1. QC -------
seurat_mtx[["percent.mt"]] <- PercentageFeatureSet(seurat_mtx, pattern = "^mt-")
seurat_mtx <- subset(seurat_mtx, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & 
                       percent.mt < 20 & nCount_RNA > 1500 & nCount_RNA < 50000)
# 17858 features across 25566 samples within 1 assay (for nfeature >200)
# 17858 features across 25555 samples (for nfeature >500)

seurat_mtx <- NormalizeData(seurat_mtx)
seurat_mtx <- FindVariableFeatures(seurat_mtx, selection.method = "vst")
all.genes <- rownames(seurat_mtx)
seurat_mtx <- ScaleData(seurat_mtx, features = all.genes)
seurat_mtx <- RunPCA(seurat_mtx, features = VariableFeatures(object = seurat_mtx))
DimHeatmap(seurat_mtx, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(seurat_mtx)
seurat_mtx <- FindNeighbors(seurat_mtx, dims = 1:13)
seurat_mtx <- FindClusters(seurat_mtx, resolution = 0.05) 
DimPlot(seurat_mtx, group.by = "RNA_snn_res.0.05", label = TRUE)
Idents(seurat_mtx) <- "RNA_snn_res.0.05"

#tSNE
seurat_mtx <- RunTSNE(seurat_mtx, dims = 1:13)
TSNEPlot(seurat_mtx, raster = F, cols =paletteer_c("grDevices::Sunset", 5))

#Option 2: UMAP
seurat_mtx <- RunUMAP(seurat_mtx, dims = 1:13)
DimPlot(seurat_mtx, reduction = "umap", label = TRUE)

# ---------------------- characterize clusters -------------------
## get cluster markers
cluster1.markers <- FindAllMarkers(seurat_mtx, logfc.threshold = .25, min.pct= 0.3, only.pos = T)
cluster1.markers <- cluster1.markers[which(cluster1.markers$p_val_adj < 0.05), ]
gl<- c()
for(i in as.character(unique(cluster1.markers$cluster))){
  temp<- cluster1.markers[which(cluster1.markers$cluster == i), ]
  temp<- temp[order(-temp$avg_log2FC), c(7, 6, 2:5)]
  cat('cluster ', i, ': ', temp$gene[1:5], '\n\n')
  gl<- c(gl, temp$gene[1:5])
  write.csv(temp, paste0("cmarks_cluster2_", i, ".csv"), row.names = F)
}

#use canonical markers to identify clusters
DotPlot(seurat_mtx, features = c("Csf1r","Cd14","Cd68","Ptprc","Tmem119","Sall1","P2ry12","Ms4a7","Adgre1","Ms4a6c","Ly6c2","Ly6g","S100a9","H2-Ab1","Cd74","Itgax","Cd3e","Ms4a4b","Klrb1c","Ly6d","Igkc","Cd19","Cldn5","Col4a1","Pdgfrb","Rgs5"),
        cols =  c("#704D9E", "#EB7999")) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text (size = 20), text = element_text(size = 20))
ggsave("dotplot_named_salvador_colors_size.png", width = 6, height = 10)

#add proposed annotation
new.clusters <- c('Microglia','Macrophages','Neutrophils','Lymphocytes','Stromal/ECs')
names(new.clusters) <- levels(seurat_mtx)
seurat_mtx <- RenameIdents(seurat_mtx, new.clusters)

#save new identity and visualize
seurat_mtx$celltype <- seurat_mtx@active.ident
TSNEPlot(seurat_mtx, raster = F, cols =paletteer_c("grDevices::Sunset", 5))
save(seurat_mtx, file = "seurat_mtx_gse205037.RData")

## ------- re-cluster myeloid -----------------------
myel <- subset(seurat_mtx, ident = "Macrophages")
myel <- ScaleData(myel, features = row.names(myel),
                  vars.to.regress = c("nFeature_RNA","nCount_RNA","percent_mito"))
myel <- FindVariableFeatures(myel, selection.method = "vst")
myel <- RunPCA(myel, npcs = 50, features = VariableFeatures(myel))
ElbowPlot(myel, 50)

myel <- FindNeighbors(myel, dims = 1:40, reduction = "pca")
myel <- FindClusters(myel, resolution = 0.4)
DimPlot(myel, group.by = "RNA_snn_res.0.4", label = TRUE)
myel <- RunTSNE(myel, dims = 1:40, verbose = T)

TSNEPlot(myel, cols =paletteer_c("grDevices::Sunset", 10))

## characterize clusters
clustmarks <- FindAllMarkers(myel, logfc.threshold = .25, min.pct= 0.3, only.pos = T)
clustmarks <- clustmarks[which(clustmarks$p_val_adj < 0.05), ]

clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10
DoHeatmap(myel, features = top10$gene,  group.colors = paletteer_c("grDevices::Sunset", 10))+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
write.csv(top10, file = "top10_salvador_2.csv")

#top10 by logfc
clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10_2
DoHeatmap(myel, features = top10_2$gene, group.colors = paletteer_c("grDevices::Sunset", 10))+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
write.csv(top10_2, file = "top10_logfc_salvador_2.csv")

top10_milich <- c("Ifi30","Gm49339","ENSMUSG00000089672","AI662270","Atp5o","Gm49342","Rab5if","Metrnl","Tmem160","Bag1","Vps28","Gm1673","Gpx3","Rtcb","Gyg","Igf1","Syngr1","Gpnmb","Pld3","Selenop","Lilr4b","Pgk1","Vps28","Pnp","Fabp5","Ftl1-ps1","Ftl1","Hmox1","Mbp","Spp1","Gm49339","F10","ENSMUSG00000089672","Cxcl3","Ccl9","Arg1","Ccr2","Cxcl2","Ltb4r1","Clec4n","Rpl23a-ps3","Gm9843","Rpl6l","Rps12-ps3","Zfas1","Timp2","Selenop","Eif3f","Rps27rt","Hexb","Pgk1","S100a6","S100a10","Eno1","Anxa2","S100a11","Lilr4b","Pfn1","Fcgr2b","Cxcl3","Pclaf","Cks1b","Birc5","Stmn1","Cdca3","Smc2","Top2a","Cdca8","Smc4","Cdk1","Lcn2","Wfdc21","Ifitm6","Serpinb1a","Pglyrp1","Pclaf","Lbr","Ngp","Camp","Trem3","Mrc1","Cd163","Lyve1","Mgl2","C4b","Fcgrt","Igfbp4","Cbr2","Cfh","Folr2","Atp6v0c","Ndufb1-ps","Lilrb4a","Rnasek","Pabpn1","Mir22hg","Gng10","Gas5","Cct6a","Timm23","Atp6v0c","Gas5","Lilrb4a","Ndufb1-ps","Rnasek","2410006H16Rik","Cct6a","Cebpd","G530011O06Rik","Snhg8")
DoHeatmap(myel, features = top10_milich)+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
save(myel, file = "myel.RData")


###reviewer's suggestions ---------
setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/")
load("myel.RData")

#Markers plots
VlnPlot(myel, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl","Sparc","Cd81","Gpr34","Serpine2","Gpnmb","Cdo1","Hpse","Ccr2","Plac8"), pt.size = 0) 
DotPlot(myel, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl","Sparc","Cd81","Gpr34","Serpine2","Gpnmb","Cdo1","Hpse","Ccr2","Plac8")) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text (size = 20), text = element_text(size = 20))
ggsave("salvador_dotplot_markers.png", width = 6, height = 10)
FeaturePlot(myel, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl","Sparc","Cd81","Gpr34","Serpine2","Gpnmb","Cdo1","Hpse","Ccr2","Plac8"))

#number and proportions of cells per cluster
cluster_counts <- table(myel@active.ident)
cluster_proportions <- prop.table(cluster_counts)
cluster_summary_salvador <- data.frame(
  Cluster = names(cluster_counts),
  Cell_Count = as.numeric(cluster_counts),
  Proportion = as.numeric(cluster_proportions)
)
save(cluster_summary_salvador, file = "cluster_summary_salvador.RData")
write.csv(cluster_summary_salvador, file = "cluster_summary_salvador.csv")

ggplot(cluster_summary, aes(x = Cluster, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Cells in Each Cluster", y = "Proportion", x = "Cluster") +
  theme_minimal()




