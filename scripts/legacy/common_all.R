library(ggvenn)

setwd("~/Google Drive/My Drive/Gensel lab/R studies/common_all/")

#top10 by adj p value
salvador_p <- c("Plin2","Lgals1","Lgals3","Fabp5","Spp1","Anxa2","Cd36","Fabp4","Lpl","Anxa1","2410006H16Rik","Rapsn","Sparc","Siglech","Cd81","Gpr34","Plxdc2","Egr1","Lag3","Serpine2","Sall1","P2ry12","Pclaf","Top2a","Hist1h1b","Mki67","Smc2","Cdk1","Tk1","Rrm1","Rrm2","Stmn1","Gpnmb","Hpse","Txn1","Pld3","Rilpl2","Plpp3","Gstm1","Ddhd1","Ahnak2","Cdo1","Birc5","Ccnb2","Cenpa","Cenpf","Ube2c","Hmmr","Cdkn3","Cdca3","Cdc20","Tpx2","H2-Aa","H2-Eb1","Plbd1","Ccr2","H2-Ab1","Plac8","Rnase6","Ciita","Klrd1","Ifitm6","Apoc4","C4b","Timp2","Ctla2a","Colec12","Lyz1","Pmepa1","Plaat3","Trf","Mrc1","Cd163","Mgl2","H2-Eb1","H2-Aa","Mrc1","H2-Ab1","Cd74","Ifitm3","C4b","Lifr","Iglc3","Ly6d","Gm21762","Ccr9","Cd8b1","Mctp2","Atp1b1","Smim5","Cd7","Atp2a1")
milich_p <- c("Rpl23a-ps3","Gm9843","Rpl13a","Rpl6l","Rps12-ps3","Lamp1","Gm8730","Bri3","Gm10076","Zfas1","Gpnmb","Ahnak2","Cd24a","Cdo1","C77080","Hpse","Lyst","Psmd1","Fabp4","Cd300lf","Lgals1","Capg","Calr","Lgals3","Pkm","Pgam1","Ldha","Pdia6","Rnh1","Spp1","Tmem119","Vps28","Ppfia4","Tmem176b","Pnp","Rasgef1b","Ptgs1","Tmem176a","Tmcc3","Fgd2","Sparc","Cd81","Tspan7","Cst7","Cst3","Pmp22","Serpine2","Xlr","Cadm1","Gpr34","Plac8","S100a4","Ifitm6","Ly6c2","Ccr2","Mcemp1","Ccl9","Hp","S100a6","S100a11","Cd163","Mrc1","Cbr2","Igfbp4","Pf4","Clec4a1","Ccl7","Marcksl1","Selenop","Fcgrt")
brennan_p <- c("Ms4a7","Fcrls","Sepp1","Pmepa1","Cx3cr1","Maf","Rnase4","Ang","Lpcat2","Tmem119","Gpnmb","Fth1","Pld3","Atp6v0d2","Fabp5","Creg1","S100a1","Txn1","Rilpl2","Hpse","Lgals1","Lgals3","Plin2","Spp1","Cstb","Cd36","Capg","Calr","Pgam1","Fabp5","Ccr2","S100a4","Ifitm3","S100a6","Cd52","Cd74","Plac8","Fxyd5","Ccl9","H2-Aa","Sparc","Serpine2","Gpr34","Cadm1","Cst7","Cd81","Plxdc2","Cd34","Siglech","Cst3","Cd79a","Ms4a1","Ly6d","Mzb1","Gimap4","Pou2af1","Cd79b","Bank1","Faim3","Cd2","Xcr1","Clec9a","Snx22","Amica1","Ifi205","F630111L10Rik","Rnase6","Plac8","Klrd1","Naaa","Rad51","Uhrf1","Ccne2","Mcm5","2810417H13Rik","Hells","Mcm2","Clspn","Mcm6","Chaf1a")

datasets_p <- c("salvador_p", "milich_p", "brennan_p")

for (i in datasets_p) {
  assign (i, unique(get(i)))
}


# Calculate intersections
intersect_salvador_milich <- intersect(salvador_p, milich_p)
intersect_salvador_brennan <- intersect(salvador_p, brennan_p)
intersect_milich_brennan <- intersect(milich_p, brennan_p)
intersect_all <- Reduce(intersect, list(salvador_p, milich_p, brennan_p))

# Create a list of intersections
intersections <- list(
  "Salvador_Milich" = intersect_salvador_milich,
  "Salvador_Brennan" = intersect_salvador_brennan,
  "Milich_Brennan" = intersect_milich_brennan,
  "All" = intersect_all
)


# Print intersections
print(intersections)
intersections <- stack(intersections)
write.csv(intersections, file = "common_genes.csv", row.names = F)


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
multiple_shared_genes <- shared_clusters %>%
  group_by(cluster_milich, cluster_brennan) %>%
  summarise(shared_genes = n()) %>%
  filter(shared_genes > 1)

