library(dplyr)

### Milich x Brennan ---------------------
##adjusted p-value
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/7dpi/")
milich_p <- read.csv("top10_pval_milich.csv")
milich_p <- milich_p[,c(7,8)]

setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
brennan_p <- read.csv("top10_brennan_adjpval_2.csv")
brennan_p <- brennan_p[,c(7,8)]

intersect_milich_brennan <- intersect(milich_p[,2], brennan_p[,2])
milich_common <- milich_p[milich_p$gene %in% intersect_milich_brennan, ]
brennan_common <- brennan_p[brennan_p$gene %in% intersect_milich_brennan, ]

shared_clusters <- merge(milich_common, brennan_common, by = "gene", suffixes = c("_milich", "_brennan"))
shared_clusters <- shared_clusters %>% arrange(cluster_milich)

multiple_shared_genes <- shared_clusters %>%
  group_by(cluster_milich, cluster_brennan) %>%
  summarise(shared_genes = n()) %>%
  filter(shared_genes > 1)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/common_all/")
write.csv(shared_clusters, file = "shared_BC_pvalue.csv")

#validate on A
setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/Macrophage clustering/")
salvador_p <- read.csv("top10_pval_salvador_2.csv")
salvador_p <- salvador_p[,c(7,8)]
salvador_common <- salvador_p[salvador_p$gene %in% intersect_milich_brennan, ]
shared_clusters_p_all <- merge(shared_clusters, salvador_common, by = "gene", suffixes = c("_milich", "_brennan", "_salvador"))
shared_clusters_p_all <- shared_clusters_p_all %>% arrange(cluster_milich)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/common_all/")
write.csv(shared_clusters_p_all, file = "shared_BCA_pvalue.csv")

##logFC
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/7dpi/")
milich_log <- read.csv("top10_logfc_milich.csv")
milich_log <- milich_log[,c(7,8)]

setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
brennan_log <- read.csv("top10_brennan_logfc_2.csv")
brennan_log <- brennan_log[,c(7,8)]

intersect_milich_brennan_log <- intersect(milich_log[,2], brennan_log[,2])
milich_common_log <- milich_log[milich_log$gene %in% intersect_milich_brennan_log, ]
brennan_common_log <- brennan_log[brennan_log$gene %in% intersect_milich_brennan_log, ]

shared_clusters_log <- merge(milich_common_log, brennan_common_log, by = "gene", suffixes = c("_milich", "_brennan"))
shared_clusters_log <- shared_clusters_log %>% arrange(cluster_milich)


setwd("~/Google Drive/My Drive/Gensel lab/R studies/common_all/")
write.csv(shared_clusters_log, file = "shared_BC_logfc.csv")

#validate on A
setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/Macrophage clustering/")
salvador_log <- read.csv("top10_logfc_salvador_2.csv")
salvador_log <- salvador_log[,c(7,8)]
salvador_common_log <- salvador_log[salvador_log$gene %in% intersect_milich_brennan_log, ]
shared_clusters_log_all <- merge(shared_clusters_log, salvador_common_log, by = "gene", suffixes = c("_milich", "_brennan", "_salvador"))
shared_clusters_log_all <- shared_clusters_log_all %>% arrange(cluster_milich)


setwd("~/Google Drive/My Drive/Gensel lab/R studies/common_all/")
write.csv(shared_clusters_log_all, file = "shared_BCA_logfc.csv")

