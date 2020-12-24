library(cola)
library(tidyverse)

df <- readRDS("data/pcawg_sample_tidy_info.rds")

mat_abs <- df %>%
  filter(keep) %>%
  select(sample, starts_with("Abs_Sig")) %>%
  column_to_rownames("sample") %>%
  t()
mat_abs <- adjust_matrix(mat_abs)
rownames(mat_abs)

mat_rel <- df %>%
  filter(keep) %>%
  select(sample, starts_with("Rel_Sig")) %>%
  column_to_rownames("sample") %>%
  t()
mat_rel <- adjust_matrix(mat_rel)
rownames(mat_rel)


# Select suitable parameters ----------------------------------------------

ds <- colSums(mat_abs)
boxplot(ds)

set.seed(123)
select_samps <- sample(names(sort(ds)), 500)

boxplot(ds[select_samps])

rl_abs <- run_all_consensus_partition_methods(mat_abs[, select_samps], top_n = 10, mc.cores = 8, max_k = 10)
cola_report(rl_abs, output_dir = "output/cola_report/pcawg_abs_sigs_500_sampls", mc.cores = 8)

rl_rel <- run_all_consensus_partition_methods(mat_rel[, select_samps], top_n = 10, mc.cores = 8, max_k = 10)
cola_report(rl_rel, output_dir = "output/cola_report/pcawg_rel_sigs_500_sampls", mc.cores = 8)

rl_cmb <- run_all_consensus_partition_methods(
  rbind(mat_abs[, select_samps],
        mat_rel[, select_samps]), top_n = 20, mc.cores = 8, max_k = 10)
cola_report(rl_cmb, output_dir = "output/cola_report/pcawg_cmb_sigs_500_sampls", mc.cores = 8)

rm(rl_abs, rl_rel, rl_cmb)

# Run clustering ----------------------------------------------------------

final_abs <- run_all_consensus_partition_methods(
  mat_abs,
  top_value_method = "ATC",
  partition_method = "skmeans",
  top_n = 10, mc.cores = 8, max_k = 10
)
cola_report(final_abs, output_dir = "output/cola_report/pcawg_abs_sigs_500_sampls_final", mc.cores = 8)

final_rel <- run_all_consensus_partition_methods(
  mat_rel,
  top_value_method = "ATC",
  partition_method = "skmeans",
  top_n = 10, mc.cores = 8, max_k = 10
)
cola_report(final_rel, output_dir = "output/cola_report/pcawg_rel_sigs_500_sampls_final", mc.cores = 8)

final_cmb <- run_all_consensus_partition_methods(
  rbind(mat_abs, mat_rel),
  top_value_method = "ATC",
  partition_method = "skmeans",
  top_n = 20, mc.cores = 8, max_k = 10
)
cola_report(final_cmb, output_dir = "output/cola_report/pcawg_cmb_sigs_500_sampls_final", mc.cores = 8)

save(final_abs, final_rel, final_cmb, file = "output/cola_report/final_result.RData")
# Check all cola reports and find ABS features is the best option
save(final_abs, file = "data/pcawg_cola_result.rds")
