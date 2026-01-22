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

#top10 by adj p value
top10_milich <- c("Rpl23a-ps3","Gm9843","Rpl13a","Rpl6l","Rps12-ps3","Lamp1","Gm8730","Bri3","Gm10076","Zfas1","Gpnmb","Ahnak2","Cd24a","Cdo1","C77080","Hpse","Lyst","Psmd1","Fabp4","Cd300lf","Lgals1","Capg","Calr","Lgals3","Pkm","Pgam1","Ldha","Pdia6","Rnh1","Spp1","Tmem119","Vps28","Ppfia4","Tmem176b","Pnp","Rasgef1b","Ptgs1","Tmem176a","Tmcc3","Fgd2","Sparc","Cd81","Tspan7","Cst7","Cst3","Pmp22","Serpine2","Xlr","Cadm1","Gpr34","Plac8","S100a4","Ifitm6","Ly6c2","Ccr2","Mcemp1","Ccl9","Hp","S100a6","S100a11","Cd163","Mrc1","Cbr2","Igfbp4","Pf4","Clec4a1","Ccl7","Marcksl1","Selenop","Fcgrt")
top10_salvador <- c("Plin2","Lgals1","Lgals3","Fabp5","Spp1","Anxa2","Cd36","Fabp4","Lpl","Anxa1","2410006H16Rik","Rapsn","Sparc","Siglech","Cd81","Gpr34","Plxdc2","Egr1","Lag3","Serpine2","Sall1","P2ry12","Pclaf","Top2a","Hist1h1b","Mki67","Smc2","Cdk1","Tk1","Rrm1","Rrm2","Stmn1","Gpnmb","Hpse","Txn1","Pld3","Rilpl2","Plpp3","Gstm1","Ddhd1","Ahnak2","Cdo1","Birc5","Ccnb2","Cenpa","Cenpf","Ube2c","Hmmr","Cdkn3","Cdca3","Cdc20","Tpx2","H2-Aa","H2-Eb1","Plbd1","Ccr2","H2-Ab1","Plac8","Rnase6","Ciita","Klrd1","Ifitm6","Apoc4","C4b","Timp2","Ctla2a","Colec12","Lyz1","Pmepa1","Plaat3","Trf","Mrc1","Cd163","Mgl2","H2-Eb1","H2-Aa","Mrc1","H2-Ab1","Cd74","Ifitm3","C4b","Lifr","Iglc3","Ly6d","Gm21762","Ccr9","Cd8b1","Mctp2","Atp1b1","Smim5","Cd7","Atp2a1")
#top10_salvador <- c("Plin2","Lgals1","Lgals3","Fabp5","Spp1","Anxa2","Cd36","Fabp4","Lpl","Anxa1","Sparc","Siglech","Cd81","Gpr34","Plxdc2","Lag3","Egr1","Serpine2","Sall1","P2ry12","2410006H16Rik","Epb41l4aos","Rapsn","Gpnmb","Pld3","Hpse","Creg1","Txn1","Rilpl2","Plpp3","Cdo1","Ahnak2","S100a1","Birc5","Ccnb2","Cenpa","Cenpf","Ube2c","Hmmr","Cdkn3","Cdca3","Cdc20","Ccnb1","Pclaf","Rrm2","Tk1","Hist1h1b","Clspn","Asf1b","Hist1h2ap","Top2a","Mcm5","E2f8","Pclaf","Top2a","Hist1h1b","Mki67","Cdk1","Prc1","Tk1","Pbk","Hist1h2ae","Rrm2","Ccr2","Ifitm6","Klrd1","Cd209a","Ciita","Plac8","H2-Aa","H2-Eb1","Plbd1","S100a4","Apoc4","C4b","Timp2","Ctla2a","Colec12","Lyz1","Pmepa1","Plaat3","Folr2","Adarb1","Cd163","Mgl2","H2-Eb1","H2-Aa","H2-Ab1","Cd74","Mrc1","Ifitm3","Tmem176b","Lifr","Plac8","Xcr1","Jaml","Ifi205","Itgae","Clec9a","Fcrla","Snx22","Gm43914","Gcsam")
top10_milich <- unique (top10_milich)
top10_salvador <- unique (top10_salvador)
top10_all <- c(top10_milich, top10_salvador)
common <- top10_all[duplicated(top10_all)]
write.csv(common, file = "common_adjpval_ultimate.csv")


##top10 by logfc
top10_2_milich <- c ("ENSMUSG00000029333","Zfas1","Rpl6l","Ctla2b","Rps26-ps1","Rpl13-ps3","Gm9843","Rpl23a-ps3","Rpl36-ps3","Septin2","C77080","Lyst","Plpp3","Npy","Ahnak2","Klhdc4","Tst","Cdo1","Ddhd1","Mfge8","Cd36","Lpl","Spp1","Ddx39","Kcnn4","Lilr4b","Pdia6","Zranb3","Calr","Stmn1","Pnp","Tmem119","Rasgef1b","Ppfia4","Vps28","Ptgs1","Tmem176b","Rogdi","Dennd4a","Tmcc3","Tspan7","Sparc","Cst7","Gpr84","Xlr","Apoc1","Serpine2","Bcl2a1a","Cadm1","Ctsf","Plac8","Ifitm6","Ly6c2","Hp","Mcemp1","S100a8","S100a4","Ccr2","Itgb7","Chil3","Cd163","Mrc1","Cbr2","Ccl12","Igfbp4","Ccl7","Pf4","Clec4a1","Hpgd","Marcksl1")
top10_2_salvador <- c("F7","Fabp4","Cd36","Il7r","Lgals3","Anxa1","Mmp8","Lpl","Mfge8","Spp1","Rapsn","2410006H16Rik","Sall1","Lag3","Nav3","Sparc","Siglech","Gpr84","Slc12a2","Cd34","P2ry12","Tanc2","Esco2","Hist1h1b","Hist1h2ab","Hist1h3c","Ccne2","Hist1h1a","Hist1h2ap","Rrm2","Hist1h2ae","Clspn","Gdf15","Fbxo32","Cdo1","Plpp3","Gadd45a","Gstm1","Rilpl2","Capn15","Hpse","Mcoln3","Cdc20","Cdkn3","Nek2","Cenpf","Ccnb1","Hmmr","Aspm","Ccnb2","Kif20a","Ube2c","Klrd1","Jaml","Ifitm6","Plac8","Ciita","H2-DMb2","Mcemp1","H2-Ab1","H2-Aa","H2-Eb1","C4b","Lyz1","Apoc4","Colec12","Adarb1","Apoc1","Siglec1","1700003F12Rik","Etv1","Folr2","Cd163","Mgl2","Mrc1","C4b","Lifr","Clec4a1","H2-Eb1","Cd74","H2-Aa","St3gal6","Klk1","Havcr1","Gm21762","Klk1b27","Iglc3","Cd8b1","Ly6d","Mctp2","Ptprf","Dntt")
#top10_2_salvador <- c("F7","Fabp4","Cd36","Il7r","Lgals3","Anxa1","Mmp8","Spp1","Lpl","Mfge8","Siglech","Sall1","Lag3","Sparc","Nav3","Gpr84","Slc12a2","Cd34","P2ry12","Tanc2","Epb41l4aos","Rapsn","2410006H16Rik","Gdf15","Fbxo32","Rragd","Cdo1","Plpp3","Capn15","Hpse","Rilpl2","Gadd45a","Gstm1","Cdc20","Cdkn3","Nek2","Cenpf","Ccnb1","Ccnb2","Hmmr","Aspm","Kif20a","Plk1","Hist1h1a","Ccne2","Hist1h1b","Ccne1","Rrm2","Hist1h2ae","Hist1h3c","Mcm10","Hist1h3f","E2f8","Esco2","Hist1h2ab","Tcf19","Hist1h3c","Hist1h1b","Neil3","Pbk","Bik","Hist1h2ap","Rad51ap1","Cd209a","Ifitm6","Klrd1","Mcemp1","Plac8","Il1b","H2-Ab1","H2-Aa","H2-Eb1","Itgb7","Lyz1","C4b","Apoc4","Colec12","Adarb1","Apoc1","1700003F12Rik","Etv1","Siglec1","Folr2","Cd163","Mgl2","Mrc1","Lifr","Clec4a1","H2-Eb1","Cd74","H2-Aa","H2-Ab1","St3gal6","Xcr1","Gcsam","Snx22","Clec9a","Ifi205","Gm43914","Mycl","Itgae","Fcrla","Cldn1")
top10_2_milich <- unique (top10_2_milich)
top10_2_salvador <- unique (top10_2_salvador)
top10_2_all <- c(top10_2_milich, top10_2_salvador)
common_2 <- top10_2_all[duplicated(top10_2_all)]
write.csv(common_2, file = "common_logfc_ultimate.csv")



## compare to brennan
#pvalue
top10_brennan <- c("Ms4a7","Fcrls","Sepp1","Pmepa1","Cx3cr1","Maf","Rnase4","Ang","Lpcat2","Tmem119","Gpnmb","Fth1","Pld3","Atp6v0d2","Fabp5","Creg1","S100a1","Txn1","Rilpl2","Hpse","Lgals1","Lgals3","Plin2","Spp1","Cstb","Cd36","Capg","Calr","Pgam1","Fabp5","Ccr2","S100a4","Ifitm3","S100a6","Cd52","Cd74","Plac8","Fxyd5","Ccl9","H2-Aa","Sparc","Serpine2","Gpr34","Cadm1","Cst7","Cd81","Plxdc2","Cd34","Siglech","Cst3","Cd79a","Ms4a1","Ly6d","Mzb1","Gimap4","Pou2af1","Cd79b","Bank1","Faim3","Cd2","Xcr1","Clec9a","Snx22","Amica1","Ifi205","F630111L10Rik","Rnase6","Plac8","Klrd1","Naaa","Rad51","Uhrf1","Ccne2","Mcm5","2810417H13Rik","Hells","Mcm2","Clspn","Mcm6","Chaf1a")
#top10_brennan <- c("Fcrls","Sepp1","Cx3cr1","Pmepa1","Aif1","Rnase4","Maf","Ang","Tmem119","Lpcat2","Gpnmb","Pld3","Fth1","Fabp5","Atp6v0d2","Creg1","S100a1","Txn1","Sgk1","Cdo1","Lgals1","Lgals3","Plin2","Cstb","Spp1","Capg","Cd36","Calr","Aprt","Pgam1","Ccr2","S100a4","Ifitm3","S100a6","Cd52","Cd74","Plac8","Ccl9","H2-Aa","Ifitm2","Sparc","Serpine2","Cadm1","Gpr34","Cst7","Cd81","Plxdc2","Cd34","Hexb","Cst3","Cd79a","Ms4a1","Mzb1","Ly6d","Gimap4","Pou2af1","Faim3","Cd79b","Bank1","Cd2","Xcr1","Clec9a","Snx22","Amica1","Ifi205","F630111L10Rik","Rnase6","Plac8","Klrd1","Naaa","Rad51","Uhrf1","Mcm5","2810417H13Rik","Mcm2","Hells","Mcm6","Chaf1a","Clspn","Gins2")
top10_brennan <- unique (top10_brennan)
all <- c(common, top10_brennan)
common_b <- all[duplicated(all)]
write.csv(common_b, file = "common_brennan_adjpvalue.csv")

#logfc
top10_brennan_2 <- c("Fcrls","Pmepa1","Tmem119","Maf","Ang","Rapsn","Rnase4","Cx3cr1","Lpcat2","Sepp1","Mfge8","S100a1","Atp6v0d2","Ppap2b","Cdo1","Gdf15","C77080","Tst","Fbxo32","Rilpl2","Cd36","F7","Ccl2","Hk3","Spp1","Plin2","Lpl","Msr1","Lgals3","Emp1","Ccr2","S100a4","Ly6c2","Ms4a4c","H2-Eb1","H2-Ab1","Plac8","H2-Aa","Ccl9","Cd74","Sparc","Cst7","Serpine2","Cadm1","Cd34","Gpr34","Ctsf","Plxdc2","Slc12a2","Cd81","Zbtb32","Ebf1","Mzb1","Pou2af1","Gimap4","Ms4a1","Blk","Cd2","Cd79a","Faim3","Xcr1","Clec9a","Snx22","F630111L10Rik","Siglech","Ifi205","Amica1","Naaa","Rnase6","Gpr171","Rad51","Ccne2","Uhrf1","2810417H13Rik","Clspn","Gins2","Apitd1","Hells","Chaf1a","Mcm5")
#top10_brennan_2 <- c("Fcrls","Pmepa1","Tmem119","Rapsn","Maf","Cx3cr1","Ang","Rnase4","Lpcat2","Sepp1","Mfge8","Atp6v0d2","S100a1","Ppap2b","C77080","Cdo1","Tst","Epas1","Fbxo32","Rnf128","Cd36","Ccl2","F7","Lpl","Plin2","Hk3","Spp1","Lgals3","Adam8","Mmp8","Ccr2","S100a4","Ly6c2","Ms4a4c","H2-Eb1","H2-Ab1","Plac8","H2-Aa","Ccl9","Ifitm3","Sparc","Serpine2","Cst7","Cadm1","Cd34","Gpr34","Plxdc2","Ctsf","P2ry12","Cd81","Zbtb32","Ebf1","Mzb1","Ms4a1","Gimap4","Pou2af1","Faim3","Cd2","Cd79a","Blk","Xcr1","Clec9a","Snx22","F630111L10Rik","Siglech","Ifi205","Amica1","Naaa","Rnase6","Gpr171","Rad51","Uhrf1","Clspn","2810417H13Rik","Hells","Chaf1a","Apitd1","Gins2","Mcm2","Cdca7")
top10_brennan_2 <- unique (top10_brennan_2)
all_2 <- c(common_2, top10_brennan_2)
common_b_2 <- all_2[duplicated(all_2)]
