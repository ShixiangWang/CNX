library(sigminer)

setwd("~/projects/CNX/")

tally_X <- readRDS("data/TCGA/tcga_cn_tally_X.rds")
tally_X_noLOH <- readRDS("data/TCGA/tcga_cn_tally_X_noLOH.rds")

# library(NMF)
# system.time(test <- NMF::nmf(t(tally_X$nmf_matrix), rank = 2, nrun = 1, .options = "tv4"))

tcga_types <- readRDS("data/TCGA/TCGA_type_info.rds")

sc <- tcga_types$cancer_type
names(sc) <- tcga_types$sample
length(unique(tcga_types$sample))

tcga_solutions <- readRDS("data/TCGA/tcga_cn_sigs_CN176_BP.rds")

tcga_expo <- bp_attribute_activity(
  tcga_solutions$object$K11,
  sample_class = sc,
  nmf_matrix = tally_X$nmf_matrix,
  return_class = "data.table",
  use_parallel = 12
)

saveRDS(tcga_expo, file = "data/TCGA/tcga_cn_sigs_CN176_activity.rds")

# hist(tcga_expo$similarity, breaks = 100)

# mat <- get_sig_similarity(pcawg_sigs, tcga_solutions$object$K11)
# pheatmap::pheatmap(mat$similarity, cluster_cols = F, cluster_rows = F)
# show_sig_profile(tcga_solutions$object$K11, style = "cosmic", mode = "copynumber", method = "X")
