library(sigminer)
library(pheatmap)
# Compare 176-components signatures ------------------------------------------

## SP
SP_TCGA <- sigprofiler_import("SP/TCGA_CN176X/", order_by_expo = TRUE)
SP_PCAWG <- sigprofiler_import("SP/PCAWG_CN176X/", order_by_expo = TRUE)

sim_SP <- get_sig_similarity(SP_PCAWG$solution, SP_TCGA$solution)
pheatmap::pheatmap(sim_SP$similarity, cluster_rows = F, cluster_cols = F, display_numbers = TRUE)

SP_TCGA2 <- sigprofiler_import("SP/TCGA_CN176X/", order_by_expo = TRUE, type = "all")
SP_PCAWG2 <- sigprofiler_import("SP/PCAWG_CN176X/", order_by_expo = TRUE, type = "all")

sim_SP <- get_sig_similarity(SP_PCAWG2$solution_list$S19, SP_TCGA$solution)
pheatmap::pheatmap(sim_SP$similarity, cluster_rows = F, cluster_cols = F, display_numbers = TRUE)
## SA
SA_TCGA <- readRDS("data/TCGA/tcga_cn_sigs_CN176_SA.rds")
SA_PCAWG <- readRDS("data/pcawg_cn_sigs_CN176_SA.rds")


SA_TCGA$K
SA_PCAWG$K

sim_SA <- get_sig_similarity(SA_PCAWG, SA_TCGA)
pheatmap::pheatmap(sim_SA$similarity, cluster_rows = F, cluster_cols = F, display_numbers = TRUE)

sim_cross <- get_sig_similarity(SP_TCGA2$solution_list$S22, SA_TCGA)
pheatmap::pheatmap(sim_cross$similarity, cluster_rows = F, cluster_cols = F, display_numbers = TRUE)

sim_cross2 <- get_sig_similarity(SP_PCAWG2$solution_list$S17, SA_PCAWG)
pheatmap::pheatmap(sim_cross2$similarity, cluster_rows = F, cluster_cols = F, display_numbers = TRUE)


# SP and BP ---------------------------------------------------------------

SP_TCGA2 <- sigprofiler_import("SP/TCGA_CN176X/", order_by_expo = TRUE, type = "all")
BP_PCAWG <- readRDS("data/pcawg_cn_sigs_CN176_signature.rds")

sig10 <- readRDS("data/pcawg_cn_sigs_CN176_10sigs.rds")
sig11 <- readRDS("data/pcawg_cn_sigs_CN176_11sigs.rds")

get_sig_similarity(sig11, sig10)
get_sig_similarity(sig11, BP_PCAWG)

sim_list <- lapply(SP_TCGA2$solution_list[-1], function(x) {
  diag(get_sig_similarity(BP_PCAWG, x)$similarity[1:11, 1:11])
})

sort(sapply(sim_list, mean), decreasing = TRUE)

sim <- get_sig_similarity(BP_PCAWG, SP_TCGA2$solution_list$S14)
pheatmap(sim$similarity, display_numbers = TRUE, cluster_rows = FALSE)


sim <- get_sig_similarity(SP_PCAWG2$solution_list$S10, SP_TCGA2$solution_list$S11)
pheatmap(sim$similarity, display_numbers = TRUE, cluster_rows = FALSE)

show_sig_profile(SP_PCAWG2$solution_list$S10, mode = "copynumber", method = "X", style = "cosmic")
show_sig_profile(SP_TCGA2$solution_list$S11, mode = "copynumber", method = "X", style = "cosmic")

sim <- get_sig_similarity(SP_PCAWG2$solution_list$S10, sig10)
pheatmap(sim$similarity, display_numbers = TRUE, cluster_rows = FALSE)

sim_list <- lapply(SP_TCGA2$solution_list, function(x) {
  diag(get_sig_similarity(SP_PCAWG2$solution_list$S10, x)$similarity[1:10, 1:10])
})

sort(sapply(sim_list, mean), decreasing = TRUE)


# sim_list2 <- lapply(SP_TCGA2$solution_list, function(x) {
#   get_sig_similarity(SP_PCAWG2$solution_list$S10, x)$similarity[1:10, 1:10]
# })
#
# sim_list2$S10

plot(10:30, sapply(sim_list, function(x) x[8]))
#
sim <- get_sig_similarity(SP_PCAWG2$solution_list$S10, BP_PCAWG)
pheatmap(sim$similarity, display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)

sim <- get_sig_similarity(SP_PCAWG2$solution_list$S14, BP_PCAWG)
pheatmap(sim$similarity, display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)

sim <- get_sig_similarity(SP_PCAWG2$solution_list$S10, SP_TCGA2$solution_list$S11)
pheatmap(sim$similarity, display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)

mean(diag(sim$similarity))

sim <- get_sig_similarity(SP_PCAWG2$solution_list$S10, SP_TCGA2$solution_list$S18)
pheatmap(sim$similarity, display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)

mean(diag(sim$similarity))

sim <- get_sig_similarity(SP_PCAWG2$solution_list$S13, SP_TCGA2$solution_list$S20)
pheatmap(sim$similarity, display_numbers = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)
mean(diag(sim$similarity))

