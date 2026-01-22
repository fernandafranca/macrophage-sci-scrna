library(patchwork)
library(Seurat)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/")
load("myel.RData") 
salvador <- myel
rm(myel)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/")
load("macro_7dpi.RData")
milich <- macro
rm(macro)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
load ("macro_brennan.RData")
brennan <- macro
rm(macro)

gene_sets <- list(
  "cluster 1" = c("Lgals1","Lgals3","Spp1","Cd36","Lpl"), 
  "cluster 2" = c("Cd81", "Gpr34", "Serpine2", "Sparc"),
  "cluster 3" = c("Gpnmb", "Hpse", "Cdo1"),
  "cluster 4" = c("Ccr2", "Plac8")
)


#Salvador----------
salvador_colors <-  c("#A1C9F4", "#BAB0AC", "#F4A3C0", "#BAB0AC", "#FFB482", "#BAB0AC", "#8DE5A1", "#BAB0AC", "#BAB0AC", "#BAB0AC")
#cluster blue
VlnPlot(salvador, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl"), cols = salvador_colors, pt.size = 0, ncol = 5) 
setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("salvador_vlnplot_bluecluster_1row.png", width = 15, height = 3)

#cluster pink
VlnPlot(salvador, features = c("Cd81", "Gpr34", "Serpine2", "Sparc"), cols = salvador_colors, pt.size = 0, ncol = 4) 
setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("salvador_vlnplot_pinkcluster.png", width = 10, height = 3)

#cluster orange
VlnPlot(salvador, features = c("Gpnmb", "Hpse", "Cdo1"), cols = salvador_colors, pt.size = 0, ncol = 3) 
ggsave("salvador_vlnplot_orangecluster.png", width = 8, height = 3)

#cluster green
VlnPlot(salvador, features = c("Ccr2", "Plac8"), cols = salvador_colors, pt.size = 0, ncol = 2) 
ggsave("salvador_vlnplot_greencluster.png", width = 5, height = 3)



#Milich
milich_colors <- c("#BAB0AC", "#FFB482", "#A1C9F4", "#BAB0AC", "#F4A3C0","#8DE5A1", "#BAB0AC")
#cluster blue
VlnPlot(milich, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl"), cols = milich_colors, pt.size = 0, ncol = 5)  
ggsave("milich_vlnplot_bluecluster_1row.png", width = 15, height = 3)
#cluster pink
VlnPlot(milich, features = c("Cd81", "Gpr34", "Serpine2", "Sparc"), cols = milich_colors, pt.size = 0, ncol = 4) 
ggsave("milich_vlnplot_pinkcluster.png", width = 10, height = 3)

#cluster orange
VlnPlot(milich, features = c("Gpnmb", "Hpse", "Cdo1"), cols = milich_colors, pt.size = 0, ncol = 3) 
ggsave("milich_vlnplot_orangecluster.png", width = 8, height = 3)

#cluster green
VlnPlot(milich, features = c("Ccr2", "Plac8"), cols = milich_colors, pt.size = 0, ncol = 2) 
ggsave("milich_vlnplot_greencluster.png", width = 5, height = 3)


#Brennan
brennan_colors <- c("#BAB0AC", "#FFB482", "#A1C9F4","#8DE5A1", "#F4A3C0", "#BAB0AC", "#BAB0AC", "#BAB0AC")
#cluster blue
VlnPlot(brennan, features = c("Lgals1","Lgals3","Spp1","Cd36","Lpl"), cols = brennan_colors, pt.size = 0, ncol = 5)  
ggsave("brennan_vlnplot_bluecluster_1row.png", width = 15, height = 3)
#cluster pink
VlnPlot(brennan, features = c("Cd81", "Gpr34", "Serpine2", "Sparc"), cols = brennan_colors, pt.size = 0, ncol = 4) 
ggsave("brennan_vlnplot_pinkcluster.png", width = 10, height = 3)

#cluster orange
VlnPlot(brennan, features = c("Gpnmb", "Hpse", "Cdo1"), cols = brennan_colors, pt.size = 0, ncol = 3) 
ggsave("brennan_vlnplot_orangecluster.png", width = 8, height = 3)

#cluster green
VlnPlot(brennan, features = c("Ccr2", "Plac8"), cols = brennan_colors, pt.size = 0, ncol = 2) 
ggsave("brennan_vlnplot_greencluster.png", width = 5, height = 3)




