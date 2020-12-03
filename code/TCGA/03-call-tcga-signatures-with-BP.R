library(sigminer)

tally_X <- readRDS("data/TCGA/tcga_cn_tally_X.rds")
tally_X_noLOH <- readRDS("data/TCGA/tcga_cn_tally_X_noLOH.rds")

library(NMF)
system.time(test <- NMF::nmf(t(tally_X$nmf_matrix), rank = 2, nrun = 1, .options = "tv4"))
