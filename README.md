{\rtf1\ansi\ansicpg1252\cocoartf2761
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\froman\fcharset0 Times-Bold;\f1\froman\fcharset0 Times-Roman;\f2\froman\fcharset0 Times-Italic;
\f3\fmodern\fcharset0 Courier;\f4\fnil\fcharset0 Menlo-Regular;\f5\fnil\fcharset0 LucidaGrande;
}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red109\green109\blue109;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;\cssrgb\c50196\c50196\c50196;}
{\*\listtable{\list\listtemplateid1\listhybrid{\listlevel\levelnfc23\levelnfcn23\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{disc\}}{\leveltext\leveltemplateid1\'01\uc0\u8226 ;}{\levelnumbers;}\fi-360\li720\lin720 }{\listname ;}\listid1}
{\list\listtemplateid2\listhybrid{\listlevel\levelnfc0\levelnfcn0\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{decimal\}}{\leveltext\leveltemplateid101\'01\'00;}{\levelnumbers\'01;}\fi-360\li720\lin720 }{\listlevel\levelnfc23\levelnfcn23\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{circle\}}{\leveltext\leveltemplateid102\'01\uc0\u9702 ;}{\levelnumbers;}\fi-360\li1440\lin1440 }{\listname ;}\listid2}
{\list\listtemplateid3\listhybrid{\listlevel\levelnfc23\levelnfcn23\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{disc\}}{\leveltext\leveltemplateid201\'01\uc0\u8226 ;}{\levelnumbers;}\fi-360\li720\lin720 }{\listname ;}\listid3}
{\list\listtemplateid4\listhybrid{\listlevel\levelnfc23\levelnfcn23\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{disc\}}{\leveltext\leveltemplateid301\'01\uc0\u8226 ;}{\levelnumbers;}\fi-360\li720\lin720 }{\listname ;}\listid4}
{\list\listtemplateid5\listhybrid{\listlevel\levelnfc23\levelnfcn23\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{disc\}}{\leveltext\leveltemplateid401\'01\uc0\u8226 ;}{\levelnumbers;}\fi-360\li720\lin720 }{\listname ;}\listid5}}
{\*\listoverridetable{\listoverride\listid1\listoverridecount0\ls1}{\listoverride\listid2\listoverridecount0\ls2}{\listoverride\listid3\listoverridecount0\ls3}{\listoverride\listid4\listoverridecount0\ls4}{\listoverride\listid5\listoverridecount0\ls5}}
\margl1440\margr1440\vieww27880\viewh18260\viewkind0
\deftab720
\pard\pardeftab720\sa240\partightenfactor0

\f0\b\fs24 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Reproducible analysis from Stapenhorst Fran\'e7a et al., Experimental Neurology (2025)
\f1\b0 \
\pard\pardeftab720\sa298\partightenfactor0

\f0\b\fs36 \cf0 Overview\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 This repository contains the analysis code and key outputs for the paper:\
\pard\pardeftab720\sa240\partightenfactor0

\f0\b \cf0 Stapenhorst Fran\'e7a F. et al. (2025)
\f1\b0 \uc0\u8232 
\f2\i Conserved macrophage subpopulations across single-cell RNA-seq datasets after spinal cord injury
\f1\i0 \uc0\u8232 
\f0\b Experimental Neurology
\f1\b0 \
In this work, we performed a cross-dataset single-cell RNA-seq analysis of intraspinal macrophages at 7 days post-injury (7 dpi) across three independent spinal cord injury datasets to identify four conserved macrophage populations (FSF1\'96FSF4) that represent reproducible immune states after injury.\
This repository is designed as a transparent, analysis-focused record of the computational workflow used to:\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa240\partightenfactor0
\ls1\ilvl0\cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 process each dataset,\
\ls1\ilvl0\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 identify macrophage clusters,\
\ls1\ilvl0\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 perform cross-dataset validation with SingleR,\
\ls1\ilvl0\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 derive top-10 gene signatures,\
\ls1\ilvl0\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 and quantify FSF population proportions.\
\pard\pardeftab720\sa240\partightenfactor0
\cf0 Raw data and large intermediate objects (e.g., Seurat objects) are intentionally excluded to keep the repository lightweight and public-safe.\
\pard\pardeftab720\partightenfactor0
\cf3 \strokec3 \
\pard\pardeftab720\sa298\partightenfactor0

\f0\b\fs36 \cf0 \strokec2 Datasets\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 We analyzed three publicly available scRNA-seq datasets at 7 dpi:\

\itap1\trowd \taflags0 \trgaph108\trleft-108 \trbrdrt\brdrnil \trbrdrl\brdrnil \trbrdrr\brdrnil 
\clvertalc \clshdrawnil \clwWidth773\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx2880
\clvertalc \clshdrawnil \clwWidth1173\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx5760
\clvertalc \clshdrawnil \clwWidth4968\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx8640
\pard\intbl\itap1\pardeftab720\qc\partightenfactor0

\f0\b \cf0 Dataset\cell 
\pard\intbl\itap1\pardeftab720\qc\partightenfactor0
\cf0 Accession\cell 
\pard\intbl\itap1\pardeftab720\qc\partightenfactor0
\cf0 Source\cell \row

\itap1\trowd \taflags0 \trgaph108\trleft-108 \trbrdrl\brdrnil \trbrdrr\brdrnil 
\clvertalc \clshdrawnil \clwWidth773\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx2880
\clvertalc \clshdrawnil \clwWidth1173\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx5760
\clvertalc \clshdrawnil \clwWidth4968\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx8640
\pard\intbl\itap1\pardeftab720\partightenfactor0
\cf0 A
\f1\b0 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0
\cf0 GSE205037\cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0
\cf0 Salvador et al.\cell \row

\itap1\trowd \taflags0 \trgaph108\trleft-108 \trbrdrl\brdrnil \trbrdrr\brdrnil 
\clvertalc \clshdrawnil \clwWidth773\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx2880
\clvertalc \clshdrawnil \clwWidth1173\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx5760
\clvertalc \clshdrawnil \clwWidth4968\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx8640
\pard\intbl\itap1\pardeftab720\partightenfactor0

\f0\b \cf0 B
\f1\b0 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0
\cf0 GSE162610\cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0
\cf0 Milich et al. (myeloid subset from authors\'92 GitHub)\cell \row

\itap1\trowd \taflags0 \trgaph108\trleft-108 \trbrdrl\brdrnil \trbrdrt\brdrnil \trbrdrr\brdrnil 
\clvertalc \clshdrawnil \clwWidth773\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx2880
\clvertalc \clshdrawnil \clwWidth1173\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx5760
\clvertalc \clshdrawnil \clwWidth4968\clftsWidth3 \clmart10 \clmarl10 \clmarb10 \clmarr10 \clbrdrt\brdrnil \clbrdrl\brdrnil \clbrdrb\brdrnil \clbrdrr\brdrnil \clpadt20 \clpadl20 \clpadb20 \clpadr20 \gaph\cellx8640
\pard\intbl\itap1\pardeftab720\partightenfactor0

\f0\b \cf0 C
\f1\b0 \cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0
\cf0 GSE196928\cell 
\pard\intbl\itap1\pardeftab720\partightenfactor0
\cf0 Brennan et al.\cell \lastrow\row
\pard\pardeftab720\sa240\partightenfactor0
\cf0 All datasets were processed using Seurat (v5) and harmonized at the macrophage level before cross-dataset comparisons.\
\pard\pardeftab720\partightenfactor0
\cf3 \strokec3 \
\pard\pardeftab720\sa298\partightenfactor0

\f0\b\fs36 \cf0 \strokec2 Scientific Goal\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 The goal of this project is to determine whether macrophage subpopulations observed after spinal cord injury are: dataset-specific artifacts\uc0\u8232 or conserved biological immune states.\
Using independent datasets and orthogonal validation (SingleR + DEG signatures), we defined four reproducible macrophage populations, or Findable Sequence Family (FSF): FSF1, FSF2. FSF3 and FSF4.\
These populations appear consistently across all three datasets.\
\pard\pardeftab720\partightenfactor0
\cf3 \strokec3 \
\pard\pardeftab720\sa298\partightenfactor0

\f0\b\fs36 \cf0 \strokec2 Repository Structure\
\pard\pardeftab720\partightenfactor0

\f3\b0\fs26 \cf0 \
macrophage-fsf-sci-scrna/\

\f4 \uc0\u9500 \u9472 \u9472 
\f3  scripts/        # All analysis scripts used for the paper\

\f4 \uc0\u9500 \u9472 \u9472 
\f3  results/        # Tables derived from the analyses\

\f4 \uc0\u9474 
\f3    
\f4 \uc0\u9500 \u9472 \u9472 
\f3  tables/     # cluster summaries, FSF gene lists, top10 DEGs\

\f4 \uc0\u9474 
\f3    
\f4 \uc0\u9492 \u9472 \u9472 
\f3  singler/    # cross-dataset SingleR comparison outputs\

\f4 \uc0\u9500 \u9472 \u9472 
\f3  figures/        # Representative figures used in the manuscript\

\f4 \uc0\u9474 
\f3    
\f4 \uc0\u9500 \u9472 \u9472 
\f3  overview/\

\f4 \uc0\u9474 
\f3    
\f4 \uc0\u9500 \u9472 \u9472 
\f3  salavdor_A/\

\f4 \uc0\u9474 
\f3    
\f4 \uc0\u9500 \u9472 \u9472 
\f3  milich_B/\

\f4 \uc0\u9474 
\f3    
\f4 \uc0\u9500 \u9472 \u9472 
\f3  brennan_C/\

\f4 \uc0\u9474 
\f3    
\f4 \uc0\u9492 \u9472 \u9472 
\f3  singler/\

\f4 \uc0\u9492 \u9472 \u9472 
\f3  README.md\
\pard\pardeftab720\sa240\partightenfactor0

\f1\fs24 \cf0 Large objects such as raw matrices, Seurat objects, and 
\f3\fs26 .RData
\f1\fs24  workspaces are excluded and ignored by 
\f3\fs26 .gitignore
\f1\fs24 .\
\pard\pardeftab720\partightenfactor0
\cf3 \strokec3 \
\pard\pardeftab720\sa298\partightenfactor0

\f0\b\fs36 \cf0 \strokec2 Analysis Workflow (high-level)\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa240\partightenfactor0
\ls2\ilvl0
\fs24 \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	1	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Process each dataset independently
\f1\b0 \
\pard\tx940\tx1440\pardeftab720\li1440\fi-1440\sa240\partightenfactor0
\ls2\ilvl1\cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	
\f5 \uc0\u9702 
\f1 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 QC, normalization, PCA, clustering\
\ls2\ilvl1\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	
\f5 \uc0\u9702 
\f1 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Identify macrophages using canonical markers\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa240\partightenfactor0
\ls2\ilvl0
\f0\b \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	2	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Macrophage subclustering
\f1\b0 \
\pard\tx940\tx1440\pardeftab720\li1440\fi-1440\sa240\partightenfactor0
\ls2\ilvl1\cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	
\f5 \uc0\u9702 
\f1 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Reclustering within each dataset to define macrophage populations\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa240\partightenfactor0
\ls2\ilvl0
\f0\b \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	3	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Cross-dataset validation
\f1\b0 \
\pard\tx940\tx1440\pardeftab720\li1440\fi-1440\sa240\partightenfactor0
\ls2\ilvl1\cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	
\f5 \uc0\u9702 
\f1 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Pairwise SingleR comparisons between datasets A/B/C\
\ls2\ilvl1\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	
\f5 \uc0\u9702 
\f1 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Statistical assessment of cluster similarity\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa240\partightenfactor0
\ls2\ilvl0
\f0\b \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	4	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Signature discovery
\f1\b0 \
\pard\tx940\tx1440\pardeftab720\li1440\fi-1440\sa240\partightenfactor0
\ls2\ilvl1\cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	
\f5 \uc0\u9702 
\f1 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Differential expression within macrophage clusters\
\ls2\ilvl1\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	
\f5 \uc0\u9702 
\f1 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Identification of top-10 gene signatures\
\ls2\ilvl1\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	
\f5 \uc0\u9702 
\f1 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Definition of FSF1\'96FSF4 gene sets\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa240\partightenfactor0
\ls2\ilvl0
\f0\b \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	5	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Population quantification
\f1\b0 \
\pard\tx940\tx1440\pardeftab720\li1440\fi-1440\sa240\partightenfactor0
\ls2\ilvl1\cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	
\f5 \uc0\u9702 
\f1 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Proportion of each FSF in each dataset\
\ls2\ilvl1\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	
\f5 \uc0\u9702 
\f1 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Statistical comparison across datasets\
\pard\pardeftab720\partightenfactor0
\cf3 \strokec3 \
\pard\pardeftab720\sa298\partightenfactor0

\f0\b\fs36 \cf0 \strokec2 Key Outputs\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 This repository includes:\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa240\partightenfactor0
\ls3\ilvl0
\f0\b \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Cluster summaries
\f1\b0  for macrophages in all datasets\
\ls3\ilvl0
\f0\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 SingleR cross-dataset similarity tables
\f1\b0 \
\ls3\ilvl0
\f0\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Top-10 gene signatures
\f1\b0  for each FSF population\
\ls3\ilvl0
\f0\b \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 Representative figures
\f1\b0  (tSNE/UMAP, dotplots, heatmaps, proportions)\
\pard\pardeftab720\sa240\partightenfactor0
\cf0 These files can be found in:\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa240\partightenfactor0
\ls4\ilvl0
\f3\fs26 \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 results/tables/
\f1\fs24 \
\ls4\ilvl0
\f3\fs26 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 results/singler/
\f1\fs24 \
\ls4\ilvl0
\f3\fs26 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 figures/
\f1\fs24 \
\pard\pardeftab720\partightenfactor0
\cf3 \strokec3 \
\pard\pardeftab720\sa298\partightenfactor0

\f0\b\fs36 \cf0 \strokec2 Reproducibility Notes\
\pard\pardeftab720\sa240\partightenfactor0

\f1\b0\fs24 \cf0 To reproduce the analysis from raw data, users should:\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa240\partightenfactor0
\ls5\ilvl0\cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 download the original datasets from GEO or authors\'92 repositories\
\ls5\ilvl0\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 use the scripts in 
\f3\fs26 scripts/
\f1\fs24 \
\ls5\ilvl0\kerning1\expnd0\expndtw0 \outl0\strokewidth0 {\listtext	\uc0\u8226 	}\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 follow the same QC and clustering parameters described in the paper\
\pard\pardeftab720\sa240\partightenfactor0
\cf0 This repository intentionally focuses on 
\f0\b analysis logic + scientific outputs
\f1\b0  rather than full data redistribution.\
\pard\pardeftab720\partightenfactor0
\cf3 \strokec3 \
\pard\pardeftab720\sa298\partightenfactor0

\f0\b\fs36 \cf0 \strokec2 Author\
\pard\pardeftab720\sa240\partightenfactor0

\fs24 \cf0 Fernanda Stapenhorst Fran\'e7a, PhD
\f1\b0 \uc0\u8232 Postdoctoral Researcher \'96 Gensel Lab\u8232 University of Kentucky\
}