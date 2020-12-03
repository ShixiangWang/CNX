library(sigminer)

tally_X <- readRDS("data/pcawg_cn_tally_X.rds")
tally_X_noLOH <- readRDS("data/pcawg_cn_tally_X_noLOH.rds")

dim(tally_X$nmf_matrix)

library(profvis)
debug(bp_extract_signatures)

sigs <- bp_extract_signatures(
  tally_X$nmf_matrix,
  range = 2:30,
  n_bootstrap = 5,
  n_nmf_run = 5,
  cores = 16,
  cores_solution = 4,
  one_batch = TRUE,
  cache_dir = "BP/BP_PCAWG_test",
  keep_cache = TRUE,
  only_core_stats = TRUE
)

bp_show_survey(sigs, fixed_ratio = F, add_score = FALSE)
saveRDS(sigs, file = "data/pcawg_cn_sigs_CN176_BP_test.rds")

library(NMF)
sigs <- sig_extract(tally_X$nmf_matrix, 10, nrun = 100, cores = 16)
saveRDS(sigs, file = "data/pcawg_cn_sigs_CN176_10sigs.rds")

sigs <- readRDS("data/pcawg_cn_sigs_CN176_10sigs.rds")
show_sig_profile(sigs, mode = "copynumber", method = "X", style = "cosmic")

sim <- get_sig_similarity(sigs, sigs)
pheatmap::pheatmap(sim$similarity, display_numbers = TRUE)
expo <- get_sig_exposure(sigs)
expo_rel <- get_sig_exposure(sigs, type = "relative")

show_cor(expo[, -1], cor_method = "pearson")
show_cor(expo_rel[, -1], cor_method = "pearson")
show_cor(expo[, -1])
show_cor(expo_rel[, -1])

plot(density(expo$Sig1))
plot(density(expo_rel$Sig1))

shapiro.test(expo$Sig1)
shapiro.test(expo_rel$Sig1)
# sigs <- bp_extract_signatures(
#   tally_X_noLOH$nmf_matrix,
#   range = 2:30,
#   n_nmf_run = 5,
#   cores = 16
# )
# saveRDS(sigs2, file = "data/pcawg_cn_sigs_CN136_BP_test.rds")
