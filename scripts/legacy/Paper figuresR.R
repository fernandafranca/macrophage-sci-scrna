library(tidyverse)
library(ggplot2)
library(Seurat)
library(ggh4x)

pastel_colors <- c("#A1C9F4", "#F4A3C0", "#FFB482", "#8DE5A1", "#D0BBFF", 
                   "#FDE051", "#FB8072", "#9DDECB", "#C6A78B", "#BAB0AC")
gradient_colors <- c("#A1C9F4", "white", "#F4A3C0")
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

## Salvador ---------
setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/")
load("seurat_mtx_gse205037.RData")

#All cells dotplot
features_salvador <-  c("Csf1r","Cd14","Cd68","Ptprc","Tmem119","Sall1","P2ry12","Ms4a7","Adgre1","Ms4a6c","Ly6c2","Ly6g","S100a9","H2-Ab1","Cd74","Itgax","Cd3e","Ms4a4b","Klrb1c","Ly6d","Igkc","Cd19","Cldn5","Col4a1","Pdgfrb","Rgs5")
DotPlot(seurat_mtx, features = features_salvador) +
        scale_color_gradient2(low = "#A1C9F4", mid = "white", high = "#F4A3C0") +
  coord_flip()+
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text (size = 15), text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))
setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("dotplot_salvador_jan25.png", width = 8, height = 10, bg = "white")

#All cells TSNE
TSNEPlot(seurat_mtx, raster = F, cols =pastel_colors)+
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
ggsave("tsne_salvador_all2.png", width = 6, height = 5, dpi = 300)

#Macrophages tSNE
setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/")
load("myel.RData")
TSNEPlot(myel, cols = pastel_colors)+
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

ggsave("tsne_salvador_macrophages2.png", width = 6, height = 5, dpi = 300)

setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")

##Heatmaps
clustmarks <- FindAllMarkers(myel, logfc.threshold = .25, min.pct= 0.3, only.pos = T)
clustmarks <- clustmarks[which(clustmarks$p_val_adj < 0.05), ]

genes_to_highlight <- unique(c("Fabp4","Lgals1","Lgals3","Spp1","Cd81","Gpr34","Serpine2","Sparc","Ahnak2","Cdo1","Gpnmb","Hpse","Ccr2","Ifitm6","Plac8","Mrc1","Cd163","Cd36","Lpl","Mfge8","Spp1","Gpr84","Sparc","Cdo1","Plpp3","Ifitm6","Mcemp1","Plac8","Apoc1","Cd163","Clec4a1","Mrc1"))

clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10

p1 <- DoHeatmap(myel, features = top10$gene,  group.colors = pastel_colors) +
  scale_fill_gradientn(colors =  rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))) 
p1 + theme(axis.text.y = element_text(color = ifelse(levels(p1$data$Feature) %in% genes_to_highlight, "red", "black")))

setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("heatmap_salvador_pvalue_2.png", width = 10.2, height = 12)

#top10 by logfc
clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10_2

p2 <- DoHeatmap(myel, features = top10_2$gene, group.colors = pastel_colors)+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
p2 + theme(axis.text.y = element_text(color = ifelse(levels(p2$data$Feature) %in% genes_to_highlight, "red", "black")))

setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("heatmap_salvador_logfc_2.png", width = 10.2, height = 12)


##Vln plot
setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/")
load("myel.RData")

#cluster1
salvador_colors <-  c("#A1C9F4", "#BAB0AC", "#F4A3C0", "#BAB0AC", "#FFB482", "#BAB0AC", "#8DE5A1", "#BAB0AC", "#BAB0AC", "#BAB0AC")
VlnPlot(myel, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl"), cols = salvador_colors, pt.size = 0) 
ggsave("salvador_vlnplot_bluecluster.png")




##Milich -------------
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/")
readRDS("myeloid-003.rds") -> myel_milich

#All cells dotplot
Idents(myel_milich) <- "sample_id"
myel_milich <- subset(myel_milich, idents = c("7dpi_sample1", "7dpi_sample2"))
Idents(myel_milich) <- "celltype"
features_milich = c('Gpr84','Ptgs1','P2ry12','Gpr34','Lag3','Siglech','Cst7', 'Ly6c2','Ccr2','F10','Plac8', 'Thbs1','Ms4a7','Pf4','Fabp4','Gpnmb', 'Top2a', 'Mki67', 'Cdk1', 'Ccnb2', "Cd74", 'S100a9','Mmp9','Ly6g','Cd177','Ltf')
DotPlot(myel_milich, features = features_milich) +
  scale_color_gradient2(low = "#A1C9F4", mid = "white", high = "#F4A3C0") +
  coord_flip() +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text (size = 15), text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))

setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("dotplot_milich.png", width = 8, height = 10, bg = "white")

#All cells TSNE
myel_milich <- RunTSNE(myel_milich, dims = 1:15, verbose = T)
TSNEPlot(myel_milich, cols = pastel_colors)+
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
ggsave("tsne_milich_all2.png", width = 6, height = 5, dpi = 300)

#Macrophages tSNE
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/")
load("macro_7dpi.RData")
TSNEPlot(macro, cols = pastel_colors)+
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
ggsave("tsne_milich_macrophages2.png", width = 6, height = 5, dpi = 300)

setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")

#Heatmaps
clustmarks <- FindAllMarkers(macro, logfc.threshold = .25, min.pct= 0.3, only.pos = T)
clustmarks <- clustmarks[which(clustmarks$p_val_adj < 0.05), ]

#top 10 by adj p value
clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10

p3 <- DoHeatmap(macro, features = top10$gene, group.colors = pastel_colors) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) )
p3 + theme(axis.text.y = element_text(color = ifelse(levels(p3$data$Feature) %in% genes_to_highlight, "red", "black")))

setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("heatmap_milich_pvalue_2.png", width = 10.2, height = 12)
#top10 by logfc
clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10_2

p4 <- DoHeatmap(macro, features = top10_2$gene, group.colors = pastel_colors) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) )
p4 + theme(axis.text.y = element_text(color = ifelse(levels(p4$data$Feature) %in% genes_to_highlight, "red", "black")))

setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("heatmap_milich_logfc_2.png", width = 10.2, height = 12)

##Brennan ------------
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
load("seurat_mtx_gse196928.RData")

#All cells dotplot
features_brennan <- c("P2ry12","Siglech","Tmem119","Gpnmb","Spp1","Fabp5","Ly6c1","Cldn5","Itm2a","H2-Aa","H2-Eb1","H2-Ab1","Slc1a2","Atp1a2","Atp1b2","Ccl5","Ms4a4b","Nkg7","Dbi","Nnat","Mt3","Ly6d","Cd79a","Cd79b","Stmn1","Hmgb2","Top2a","S100a9","S100a8","Retnlg","Hba-a2","Hbb-bs","Alas2","Olig1","Cspg5","Vcan","Apod","Dcn","Vtn","Mgp","Igfbp6","Gsn")
DotPlot(seurat_mtx, features = features_brennan) +
        scale_color_gradient2(low = "#A1C9F4", mid = "white", high = "#F4A3C0") +
          coord_flip() +
          theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1), axis.text.y = element_text (size = 15), text = element_text(size = 15)) +
          guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2)) +
          labs (color = "Average Expression", size = "Percent Expressed")
setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("dotplot_brennan_jan25.png", width = 6, height = 10, bg = "white")

#All cells TSNE
TSNEPlot(seurat_mtx, raster = F, cols =pastel_colors)+
  guides (x = axis, y = axis, color = guide_legend(nrow = 5, override.aes = list(size = 3)))+
  labs(x = "tSNE 1", y = "tSNE 2", title = "Spinal cord cells")+
  theme(axis.title = element_text(hjust = 0),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        plot.title = element_text(size = 20, face = "plain", hjust = 0.5))
ggsave("tsne_brennan_all2.png", width = 6, height = 7, dpi = 300)

#Macrophages tSNE
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
load ("macro_brennan.RData")

TSNEPlot(macro, cols = pastel_colors)+
  guides (x = axis, y = axis, color = guide_legend(nrow = 2, byrow = T, override.aes = list(size = 3)))+
  labs(x = "tSNE 1", y = "tSNE 2", title = "Macrophages") +
  theme(axis.title = element_text(hjust = 0),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.5, "cm"),
        legend.key.width = unit(1, "cm"),
        plot.title = element_text(size = 20, face = "plain", hjust = 0.5))
ggsave("tsne_brennan_macrophages2.png", width = 6, height = 5, dpi = 300)


setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")

#Heatmaps
clustmarks <- FindAllMarkers(macro, logfc.threshold = .25, min.pct= 0.3, only.pos = T)
clustmarks <- clustmarks[which(clustmarks$p_val_adj < 0.05), ]

#top 10 by adj p value
clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10

p5 <- DoHeatmap(macro, features = top10$gene,  group.colors = pastel_colors) +
  scale_fill_gradientn(colors =  rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
p5 + theme(axis.text.y = element_text(color = ifelse(levels(p5$data$Feature) %in% genes_to_highlight, "red", "black")))

setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("heatmap_brennan_pvalue_2.png", width = 10.2, height = 12)

#top 10 by logfc
clustmarks %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head (n = 10) %>%
  ungroup () -> top10_2

p6 <- DoHeatmap(macro, features = top10_2$gene,  group.colors = pastel_colors) +
  scale_fill_gradientn(colors =  rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
p6 + theme(axis.text.y = element_text(color = ifelse(levels(p6$data$Feature) %in% genes_to_highlight, "red", "black")))

setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("heatmap_brennan_logfc_2.png", width = 10.2, height = 12)
