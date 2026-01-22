🧬 Macrophage FSF Atlas in Spinal Cord Injury (scRNA-seq)

Reproducible analysis from Stapenhorst França et al., Experimental Neurology (2025)

📖 Overview

This repository contains the analysis code and key outputs for the paper:

Stapenhorst França F. et al. (2025)
Conserved macrophage subpopulations across single-cell RNA-seq datasets after spinal cord injury
Experimental Neurology

In this work, we performed a cross-dataset single-cell RNA-seq analysis of intraspinal macrophages at 7 days post-injury (7 dpi) across three independent spinal cord injury datasets to identify four conserved macrophage populations (FSF1–FSF4) that represent reproducible immune states after injury.

This repository is designed as a transparent, analysis-focused record of the computational workflow used to:

process each dataset,

identify macrophage clusters,

perform cross-dataset validation with SingleR,

derive top-10 gene signatures,

and quantify FSF population proportions.

Raw data and large intermediate objects (e.g., Seurat objects) are intentionally excluded to keep the repository lightweight and public-safe.

🧪 Datasets

We analyzed three publicly available scRNA-seq datasets at 7 dpi:

Dataset	Accession	Source
A	GSE205037	Salvador et al.
B	GSE162610	Milich et al. (myeloid subset from authors’ GitHub)
C	GSE196928	Brennan et al.

All datasets were processed using Seurat (v5) and harmonized at the macrophage level before cross-dataset comparisons.

🧬 Scientific Goal

The goal of this project is to determine whether macrophage subpopulations observed after spinal cord injury are dataset-specific artifacts or conserved biological immune states.

Using independent datasets and orthogonal validation (SingleR + DEG signatures), we defined four reproducible macrophage populations, or Findable Sequence Family (FSF): FSF1, FSF2, FSF3 and FSF4.


These populations appear consistently across all three datasets.

📁 Repository Structure
macrophage-sci-scrna/
├── scripts/        # All analysis scripts used for the paper
├── results/        # Tables derived from the analyses
│   ├── tables/     # cluster summaries, FSF gene lists, top10 DEGs
│   └── singler/    # cross-dataset SingleR comparison outputs
├── figures/        # Representative figures used in the manuscript
└── README.md


Large objects such as raw matrices, Seurat objects, and .RData workspaces are excluded and ignored by .gitignore.

🧠 Analysis Workflow (high-level)

Process each dataset independently

QC, normalization, PCA, clustering

Identify macrophages using canonical markers

Macrophage subclustering

Reclustering within each dataset to define macrophage populations

Cross-dataset validation

Pairwise SingleR comparisons between datasets A/B/C

Statistical assessment of cluster similarity

Signature discovery

Differential expression within macrophage clusters

Identification of top-10 gene signatures

Definition of FSF1–FSF4 gene sets

Population quantification

Proportion of each FSF in each dataset

Statistical comparison across datasets

📊 Key Outputs

This repository includes:

Cluster summaries for macrophages in all datasets

SingleR cross-dataset similarity tables

Top-10 gene signatures for each FSF population

Representative figures (tSNE/UMAP, dotplots, heatmaps, proportions)

These files can be found in:

results/tables/

results/singler/

figures/

🔬 Reproducibility Notes

To reproduce the analysis from raw data, users should:

download the original datasets from GEO or authors’ repositories

use the scripts in scripts/

follow the same QC and clustering parameters described in the paper

This repository intentionally focuses on analysis logic + scientific outputs rather than full data redistribution.

👩‍🔬 Author

Fernanda Stapenhorst França, PhD
Postdoctoral Researcher – Gensel Lab
University of Kentucky