library(dplyr)

### Salvador x Milich ---------------------
##adjusted p-value
setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/Macrophage clustering/")
salvador_p <- read.csv("top10_pval_salvador_2.csv")
salvador_p <- salvador_p[,c(3,6,7,8)]

setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/7dpi/")
milich_p <- read.csv("top10_pval_milich.csv")
milich_p <- milich_p[,c(3,6,7,8)]

intersect_salvador_milich <- intersect(salvador_p[,4], milich_p[,4])
salvador_common <- salvador_p[salvador_p$gene %in% intersect_salvador_milich, ]
milich_common <- milich_p[milich_p$gene %in% intersect_salvador_milich, ]

shared_clusters <- merge(salvador_common, milich_common, by = "gene", suffixes = c("_salvador", "_milich"))
shared_clusters <- shared_clusters %>% select(gene, starts_with("cluster"), starts_with("p_val_adj"), starts_with("avg_log2FC"))
shared_clusters <- shared_clusters %>% arrange(cluster_salvador)

multiple_shared_genes <- shared_clusters %>%
  group_by(cluster_salvador, cluster_milich) %>%
  summarise(shared_genes = n()) %>%
  filter(shared_genes > 1)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/common_all/")
write.csv(shared_clusters, file = "shared_AB_pvalue_2.csv")

#Validate on C
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
brennan_p <- read.csv("top10_brennan_adjpval_2.csv")
brennan_p <- brennan_p[,c(3,6,7,8)]

brennan_common <- brennan_p[brennan_p$gene %in% intersect_salvador_milich, ]
shared_clusters_p_all <- merge(shared_clusters, brennan_common, by = "gene", suffixes = c("_salvador", "_milich", "_brennan"))
shared_clusters_p_all <- shared_clusters_p_all %>% 
  select(gene, starts_with("cluster"), starts_with("p_val_adj"), starts_with("avg_log2FC")) %>%
  arrange(cluster_salvador) %>%
  mutate(across(starts_with("avg_log2FC"), ~ round(.x, 2)))
  

setwd("~/Google Drive/My Drive/Gensel lab/R studies/common_all/")
write.csv(shared_clusters_p_all, file = "shared_ABC_pvalue2.csv")

##logFC
setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/Macrophage clustering/")
salvador_log <- read.csv("top10_logfc_salvador_2.csv")
salvador_log <- salvador_log[,c(3, 6, 7,8)]

setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/7dpi/")
milich_log <- read.csv("top10_logfc_milich.csv")
milich_log <- milich_log[,c(3, 6, 7,8)]

intersect_salvador_milich_log <- intersect(salvador_log[,4], milich_log[,4])
salvador_common_log <- salvador_log[salvador_log$gene %in% intersect_salvador_milich_log, ]
milich_common_log <- milich_log[milich_log$gene %in% intersect_salvador_milich_log, ]

shared_clusters_log <- merge(salvador_common_log, milich_common_log, by = "gene", suffixes = c("_salvador", "_milich"))
shared_clusters_log <- shared_clusters_log %>% 
  select(gene, starts_with("cluster"), starts_with("p_val_adj"), starts_with("avg_log2FC")) %>%
  arrange(cluster_salvador) %>%
  mutate(across(starts_with("avg_log2FC"), ~ round(.x, 2)))

setwd("~/Google Drive/My Drive/Gensel lab/R studies/common_all/")
write.csv(shared_clusters_log, file = "shared_AB_logfc2.csv")

#Validate on C
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
brennan_log <- read.csv("top10_brennan_logfc_2.csv")
brennan_log <- brennan_log[,c(3, 6, 7,8)]

brennan_common_log <- brennan_log[brennan_log$gene %in% intersect_salvador_milich_log, ]
shared_clusters_log_all <- merge(shared_clusters_log, brennan_common_log, by = "gene", suffixes = c("_salvador", "_milich", "_brennan"))
shared_clusters_log_all <- shared_clusters_log_all %>% 
  select(gene, starts_with("cluster"), starts_with("p_val_adj"), starts_with("avg_log2FC")) %>%
  arrange(cluster_salvador) %>%
  mutate(across(starts_with("avg_log2FC"), ~ round(.x, 2)))

setwd("~/Google Drive/My Drive/Gensel lab/R studies/common_all/")
write.csv(shared_clusters_log_all, file = "shared_ABC_log2.csv")

