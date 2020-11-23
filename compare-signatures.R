library(sigminer)
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
