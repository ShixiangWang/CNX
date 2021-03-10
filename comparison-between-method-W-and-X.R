# Use PCAWG Breast and OV data as example datasets to
# compare method W and X in
# 1) reconstructed profile RSS and cosine similarity
# 2) the difference between estimated segments and observed segments
# 3) the mean similarity within copy number signatures
#
# Here we set the signature number from 3 to 10.

# Prepare data ------------------------------------------------------------

library(sigminer)
library(tidyverse)

tally_X <- readRDS("data/pcawg_cn_tally_X.rds")
pcawg_types <- readRDS("data/pcawg_type_info.rds")

cn_obj <- readRDS("data/pcawg_cn_obj.rds")
tally_W <- sig_tally(cn_obj,
                     method = "W",
                     cores = 10)
saveRDS(tally_W, file = "data/pcawg_cn_tally_W.rds")

# Extract specified signatures --------------------------------------------

ResultX <- list()
ResultW <- list()
types <- c("Breast", "Ovary-AdenoCA", "Prost-AdenoCA", "Liver-HCC", "Skin-Melanoma")

for (type in types) {
  message("=> Handling type: ", type)
  samples = pcawg_types$sample[pcawg_types$cancer_type == type]
  matX = tally_X$nmf_matrix[samples, ]
  matW = tally_W$nmf_matrix[samples, ]

  sigX = map(3:10, ~sig_extract(matX, n_sig = ., nrun = 50, cores = 10))
  sigW = map(3:10, ~sig_extract(matW, n_sig = ., nrun = 50, cores = 10))

  ResultX[[type]] = sigX
  ResultW[[type]] = sigW
  rm(sigX, sigW)
}

save(ResultX, ResultW, file = "Comparison_Sig_Data.RData")
