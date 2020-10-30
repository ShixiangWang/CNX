library(sigminer)

tally_X <- readRDS("data/pcawg_cn_tally_X.rds")
tally_X_noLOH <- readRDS("data/pcawg_cn_tally_X_noLOH.rds")

dim(tally_X$nmf_matrix)

sigs <- sig_auto_extract(
  tally_X$nmf_matrix,
  result_prefix = "CN176",
  destdir = "SA",
  K0 = 40,
  nrun = 100,
  strategy = "stable",
  cores = 20,
  skip = TRUE)
saveRDS(sigs, file = "data/pcawg_cn_sigs_CN176_SA.rds")

sigs2 <- sig_auto_extract(
  tally_X_noLOH$simplified_matrix,
  result_prefix = "CN136",
  destdir = "SA",
  K0 = 40,
  nrun = 100,
  strategy = "stable",
  cores = 16,
  skip = TRUE)
saveRDS(sigs2, file = "data/pcawg_cn_sigs_CN136_SA.rds")
