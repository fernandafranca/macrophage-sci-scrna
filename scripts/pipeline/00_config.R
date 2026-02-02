# scripts/pipeline/00_config.R
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

set.seed(123)

# repo-relative paths (assumes you run from repo root)
DIR_RESULTS <- file.path("results")
DIR_TABLES  <- file.path(DIR_RESULTS, "tables")
DIR_FIGS    <- file.path("figures")
DIR_FIG_OV  <- file.path(DIR_FIGS, "overview")

dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TABLES,  showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIGS,    showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_FIG_OV,  showWarnings = FALSE, recursive = TRUE)

# local-only data folder (NOT tracked by git)
DIR_LOCAL_DATA <- file.path("local_data")

# local-only cache for big objects (ignored by git)
DIR_LOCAL_CACHE <- file.path("local_cache")
dir.create(DIR_LOCAL_CACHE, showWarnings = FALSE, recursive = TRUE)

message("Config loaded. Working directory: ", getwd())
