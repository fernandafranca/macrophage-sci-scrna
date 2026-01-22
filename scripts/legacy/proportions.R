library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/")
load("cluster_summary_salvador.RData")
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/")
load("cluster_summary_milich.RData")
setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE196928/")
load("cluster_summary_brennan.RData")

cluster_summary_salvador <- cluster_summary_salvador %>%
  mutate(Cluster = case_when(
    Cluster == "0" ~ "FSF1",
    Cluster == "2" ~ "FSF2",
    Cluster == "4" ~ "FSF3",
    Cluster == "6" ~ "FSF4",
    TRUE ~ Cluster  # Keep other names unchanged
  ))

cluster_summary_milich <- cluster_summary_milich %>%
  mutate(Cluster = case_when(
    Cluster == "2" ~ "FSF1",
    Cluster == "4" ~ "FSF2",
    Cluster == "1" ~ "FSF3",
    Cluster == "5" ~ "FSF4",
    TRUE ~ Cluster  # Keep other names unchanged
  ))

cluster_summary_brennan <- cluster_summary_brennan %>%
  mutate(Cluster = case_when(
    Cluster == "2" ~ "FSF1",
    Cluster == "4" ~ "FSF2",
    Cluster == "1" ~ "FSF3",
    Cluster == "3" ~ "FSF4",
    TRUE ~ Cluster  # Keep other names unchanged
  ))

cluster_summary_salvador <- cluster_summary_salvador[c(1, 3, 5, 7), ]
cluster_summary_milich <- cluster_summary_milich[c(2, 3, 5, 6), ]
cluster_summary_brennan <- cluster_summary_brennan[c(2, 3, 4, 5), ]

cluster_summary_milich <- cluster_summary_milich %>% rename("Cell_Count_Milich" = "Cell_Count")
cluster_summary_brennan <- cluster_summary_brennan %>% rename("Cell_Count_Brennan" = "Cell_Count")
cluster_summary_salvador <- cluster_summary_salvador %>% rename("Cell_Count_Salvador" = "Cell_Count")

cluster_summary_milich <- cluster_summary_milich %>% rename("Proportion_Milich" = "Proportion")
cluster_summary_brennan <- cluster_summary_brennan %>% rename("Proportion_Brennan" = "Proportion")
cluster_summary_salvador <- cluster_summary_salvador %>% rename("Proportion_Salvador" = "Proportion")

combined_summary <- cbind(
  cluster_summary_salvador,
  cluster_summary_milich[, -1],
  cluster_summary_brennan[, -1]  # Remove duplicate 'Cluster' column
  
)

# Reshape proportions into long format
combined_summary_long <- pivot_longer(
  combined_summary, 
  cols = c(Proportion_Milich, Proportion_Brennan, Proportion_Salvador), 
  names_to = "Dataset", 
  values_to = "Proportion"
)

# Clean dataset names
combined_summary_long$Dataset <- gsub("Proportion_", "", combined_summary_long$Dataset)
combined_summary_long <- combined_summary_long %>% 
                      mutate (Dataset = recode(Dataset,
                                               "Salvador" = "Dataset A",
                                               "Milich" = "Dataset B",
                                               "Brennan" = "Dataset C"))

# Reshape cell counts into long format
cell_counts_long <- pivot_longer(
  combined_summary, 
  cols = c(Cell_Count_Milich, Cell_Count_Brennan, Cell_Count_Salvador), 
  names_to = "Dataset", 
  values_to = "Cell_Count"
)

# Clean dataset names
cell_counts_long$Dataset <- gsub("Cell_Count_", "", cell_counts_long$Dataset)
cell_counts_long <- cell_counts_long %>% 
  mutate (Dataset = recode(Dataset,
                           "Salvador" = "Dataset A",
                           "Milich" = "Dataset B",
                           "Brennan" = "Dataset C"))


# Merge the two long tables to ensure correct mapping
combined_summary_long <- left_join(combined_summary_long, cell_counts_long, by = c("Cluster", "Dataset"))

#colors <-  c("#A1C9F4", "#F4A3C0", "#FFB482", "#8DE5A1")

# Plot
ggplot(combined_summary_long, aes(x = Cluster, y = Proportion, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Cell_Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 4, color = "black") +
  labs(title = "Cluster Proportions Across Datasets", y = "Proportion", x = "Cluster") +
  theme_minimal()

ggplot(combined_summary_long, aes(x=Dataset, y=Proportion, fill=Cluster)) +
  geom_bar(stat = "identity") +
  labs(title = "Cluster Proportions Across Datasets", y = "Proportion", x = "Dataset") +
  scale_fill_manual(values = c("FSF1" = "#A1C9F4",
                               "FSF2" = "#F4A3C0",
                               "FSF3" = "#FFB482",
                               "FSF4" = "#8DE5A1"
    
  ))+
  theme_minimal()+
  theme (panel.grid.minor = element_blank())
setwd("~/Google Drive/My Drive/Gensel lab/macrophages review/Figures/Jan 25 Figures")
ggsave("proportions.png", width = 8, height = 12, bg = "white")
