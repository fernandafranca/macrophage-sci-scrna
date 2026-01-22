library(tidyverse)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(scran)
library(SeuratDisk)
library(HDF5Array)
library(zinbwave)
library(BiocNeighbors)
library(clusterProfiler)
library(enrichplot)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/GSE162610_Milich/")

load("macro_seurat_gse162610.RData")
load("~/Google Drive/My Drive/Gensel lab/R studies/Salvador2023/myel3.RData")

##top 10 by adjusted p value
top10_milich <- c("Gm49339","ENSMUSG00000089672","Ifi30","Atp5o","AI662270","Metrnl","Bri3","Rab5if","Erh","Vapa","Lilr4b","Pgk1","Vps28","Pnp","Mbp","Timm8a1","Chmp1b","Rasgef1b","Sbno2","Creld2","Eif3f","Rpl13a","Apoe","Rps27","Selenop","Gm9843","Rpl23a-ps3","Rpl6l","Rps27rt","Gnas","Rpl15","Apoe","Cox7a2l","Vps28","Gyg","Fosb","Rtcb","Mbp","Rpl36a-ps1","Dennd4a","Mcm5","Ccnd1","Pclaf","Stmn1","Dut","Lig1","Mcm2","Mcm6","Tipin","Gins2","Plac8","Ifitm6","Dmkn","Ly6c2","Ifitm3","Ifitm2","Vcan","S100a6","S100a11","Srgn")
top10_salvador <- c("Adam8","Clec4d","Hmox1","Nme1","Capg","Pgam1","Plin2","Pkm","Cstb","Lgals1","Apoe","Selenop","Cybb","Gyg","Ypel3","Rtcb","2410006H16Rik","Ramp1","Gm1673","AU020206","Plxdc2","Gpr34","Serpine2","Siglech","Sparc","Cd34","Cd81","Nav3","Sall1","Tspan7","H2-Eb1","Ciita","Klrd1","H2-Aa","H2-Oa","H2-DMb2","Kmo","Jaml","H2-Ab1","Ifitm1","Sell","Plac8","Ear2","Ifitm6","C3","Gpr141","Hp","Ly6c2","Vcan","Pglyrp1","Cd163","Vcam1","Cp","Lilra5","Mgl2","Bank1","Mrc1","C4b","Sult1a1","Slc15a2")
top10_milich <- unique (top10_milich)
top10_salvador <- unique (top10_salvador)
top10_all <- c(top10_milich, top10_salvador)
common <- top10_all[duplicated(top10_all)]
write.csv(common, file = "common_adjpval_3dpi.csv")

##top10 arranged from bigger logfc to smaller per cluster
top10_2_milich <- c ("ENSMUSG00000089672","Gm49339","AI662270","Gm49342","Fam177a","Ifi30","Gclm","Fdx2","Srp54c","Glrx5","Lilr4b","Pgk1","Timm8a1","Chmp1b","Pnp","Vps28","Rasgef1b","Mbp","Sbno2","Ly6c2","Zfas1","Ctla2b","Klf2","Zfp622","Eif3f","Rpl6l","Btg2","Il1b","Tmem176b","Cx3cr1","Rpl36a-ps1","Dennd4a","Rapsn","Fosb","Sbno2","Vps28","Rpl15","Gyg","Mpst","Ldhb","Mcm5","Pclaf","Ccnd1","Dut","Mcm2","Mcm6","Lig1","Stmn1","Ccl12","Gins2","Plac8","Ifitm6","Vcan","Dmkn","Hp","Klra2","Ltb4r1","Ly6c2","Mcemp1","Ms4a4c")
top10_2_salvador <- c("Ccl7","Arg1","Hmox1","Ccl2","Pdpn","Clec4d","Slc7a11","Thbs1","Rgcc","Pf4","Rapsn","Gm1673","Gpx3","Gyg","Apoe","Zfp704","Selenop","Axl","Ramp1","Crebrf","Nav3","Sall3","Sparc","Atp8a2","Slc2a5","Adamts1","Sall1","Siglech","Col27a1","Capn3","Klrd1","Ifitm1","P2ry10","Kmo","Klrb1b","Dpp4","Kctd14","H2-Oa","H2-DMb2","Il1r2","Ear2","Sell","Pglyrp1","Spn","Ly6c2","Hp","Ccdc88c","Plac8","C3","Nfe2","Vcam1","Cd163","Cp","Slc15a2","Lilra5","C4b","Mrc1","Bank1","Sult1a1","Kitl")
top10_2_milich <- unique (top10_2_milich)
top10_2_salvador <- unique (top10_2_salvador)
top10_2_all <- c(top10_2_milich, top10_2_salvador)
common_2 <- top10_2_all[duplicated(top10_2_all)]
write.csv(common_2, file = "common_logfc_3dpi.csv")
