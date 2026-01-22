library(tidyverse)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(SeuratDisk)
library(HDF5Array)
library(paletteer)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/")
readRDS("myeloid-003.rds") -> myel

colpallete <- c("#5356CA", "#A2D2FF", "#FFC8DD", "#AF9CFF", "#26547C",  "#F4BAEA",  "#FFAFCC","#E6BEDC","#CDB4DB")
#subsetting for 7 dpi
Idents(myel) <- "sample_id"
myel <- subset(myel, idents = c("7dpi_sample1", "7dpi_sample2"))
Idents(myel) <- "celltype"
ElbowPlot(myel, 50)
myel <- RunTSNE(myel, dims = 1:15, verbose = T)
TSNEPlot(myel, cols =paletteer_c("grDevices::Sunset", 6))

save(myel, file = "myel_milich.RData")
#use canonical markers to identify clusters
DotPlot(myel, features = c('Gpr84','Ptgs1','P2ry12','Gpr34','Lag3','Siglech','Cst7', 'Ly6c2','Ccr2','F10','Plac8', 'Thbs1','Ms4a7','Pf4','Fabp4','Gpnmb', 'Top2a', 'Mki67', 'Cdk1', 'Ccnb2', "Cd74", 'S100a9','Mmp9','Ly6g','Cd177','Ltf'), 
        cols =  c("#704D9E", "#EB7999")) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text (size = 20), text = element_text(size = 20))
ggsave("dotplot_named_milich_size.png", width = 6, height = 10)


#subsetting for macrophages only
macro <- subset (myel, idents = c("Macrophage", "Monocyte"))
macro <- ScaleData(macro, features = row.names(macro),
                  vars.to.regress = c("nFeature_RNA","nCount_RNA","percent_mito"))
macro <- FindVariableFeatures(macro, selection.method = "vst")
macro <- RunPCA(macro, npcs = 50, features = VariableFeatures(macro))
ElbowPlot(macro, 50)

macro <- FindNeighbors(macro, dims = 1:30, reduction = "pca")
macro <- FindClusters(macro, resolution = 0.31)
DimPlot(macro, group.by = "RNA_snn_res.0.31", label = TRUE)
macro <- RunTSNE(macro, dims = 1:40, verbose = T)

TSNEPlot(macro, cols =paletteer_c("grDevices::Sunset", 7))
DimPlot(macro, reduction = "umap", label = TRUE)

clustmarks <- FindAllMarkers(macro, logfc.threshold = .25, min.pct= 0.3, only.pos = T)
clustmarks <- clustmarks[which(clustmarks$p_val_adj < 0.05), ]

gl<- c()
for(i in as.character(unique(clustmarks$cluster))){
  temp<- clustmarks[which(clustmarks$cluster == i), ]
  temp<- temp[order(-temp$avg_log2FC), c(7, 6, 2:5)]
  cat('cluster ', i, ': ', temp$gene[1:5], '\n\n')
  gl<- c(gl, temp$gene[1:5])
  write.csv(temp, paste0("macs_cmarks_cluster_7dpi_", i, ".csv"), row.names = F)
}


##get top 10 genes
#top 10 by adj p value
clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10
DoHeatmap(macro, features = top10$gene, group.colors = paletteer_c("grDevices::Sunset", 7)) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) )

write.csv(top10, file = "top10_7dpi_withmonocytes.csv")


#top10 by logfc
clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10_2
DoHeatmap(macro, features = top10_2$gene, group.colors = paletteer_c("grDevices::Sunset", 7)) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))

write.csv(top10_2, file = "top10_7dpi_arranged_withmonocytes.csv")
save(macro, file = "macro_7dpi.RData")

###reviewer's suggestions ---------
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/")
load("macro_7dpi.RData")

#Markers plots
VlnPlot(macro, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl","Sparc","Cd81","Gpr34","Serpine2","Gpnmb","Cdo1","Hpse","Ccr2","Plac8"), pt.size = 0) 
DotPlot(macro, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl","Sparc","Cd81","Gpr34","Serpine2","Gpnmb","Cdo1","Hpse","Ccr2","Plac8")) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text (size = 20), text = element_text(size = 20))
ggsave("milich_dotplot_markers.png", width = 6, height = 10)
FeaturePlot(macro, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl","Sparc","Cd81","Gpr34","Serpine2","Gpnmb","Cdo1","Hpse","Ccr2","Plac8"), reduction = "tsne")

#number and proportions of cells per cluster
cluster_counts <- table(macro@active.ident)
cluster_proportions <- prop.table(cluster_counts)
cluster_summary_milich <- data.frame(
  Cluster = names(cluster_counts),
  Cell_Count = as.numeric(cluster_counts),
  Proportion = as.numeric(cluster_proportions)
)
save(cluster_summary_milich, file = "cluster_summary_milich.RData")
write.csv(cluster_summary_milich, file = "cluster_summary_milich.csv")

ggplot(cluster_summary_milich, aes(x = Cluster, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity") +
  labs(title = "Proportion of Cells in Each Cluster", y = "Proportion", x = "Cluster") +
  theme_minimal()



