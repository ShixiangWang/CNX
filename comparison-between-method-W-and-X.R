# Use PCAWG Breast and OV data as example datasets to
# compare method W and X in
# 1) the difference between estimated segments and observed segments
# 2) the mean similarity within copy number signatures
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


# Comparison --------------------------------------------------------------

# within similarity
get_within_similarity <- function(x) {
  sim <- suppressMessages(get_sig_similarity(x, x, set_order = FALSE))$similarity
  sim[upper.tri(sim)]
}

CrossSimX <- map_df(ResultX, function(x) {
  map(x, get_within_similarity) %>% set_names(paste0(3:10, "_Sigs")) %>%
    enframe("sig_type", "sim") %>%
    unnest("sim")
}, .id = "cancer_type")


CrossSimW <- map_df(ResultW, function(x) {
  map(x, get_within_similarity) %>% set_names(paste0(3:10, "_Sigs")) %>%
    enframe("sig_type", "sim") %>%
    unnest("sim")
}, .id = "cancer_type")

CrossSim <- bind_rows(
  CrossSimX %>% mutate(method = "X"),
  CrossSimW %>% mutate(method = "W")
)

library(ggpubr)
p <- ggboxplot(CrossSim %>%
                 mutate(
                   sig_type = stringr::str_replace(sig_type, "_Sigs", "")
                 ),
               x = "sig_type", y = "sim", fill = "method", facet.by = "cancer_type",
               xlab = "Extracted signature number", ylab = "Similarity between each other") +
  stat_compare_means(aes(group = method, label = ..p.signif..))
ggsave("output/similarity_comparison_between_methods.pdf", plot = p, width = 12, height = 7)

cn_obj <- readRDS("data/pcawg_cn_obj.rds")

observed_count <- cn_obj@summary.per.sample[, .(sample, n_of_seg)]

get_estimated_count <- function(x) {
  get_sig_exposure(x) %>%
    mutate(count = rowSums(.[, -1])) %>%
    select(sample, count)
}

count_df <- map2_df(ResultX, ResultW, function(x, y) {
  countX <- map(x, get_estimated_count) %>% set_names(paste0(3:10, "_Sigs")) %>%
    data.table::rbindlist(idcol = "sig_type")
  countX$method = "X"
  countW <- map(y, get_estimated_count) %>% set_names(paste0(3:10, "_Sigs")) %>%
    data.table::rbindlist(idcol = "sig_type")
  countW$method = "W"
  rbind(countX, countW)
}, .id = "cancer_type")

count_df2 <- left_join(count_df, observed_count, by = "sample")

save(CrossSim, count_df2, file = "Comparison_Result_Data.RData")

ggscatter(count_df2 %>% filter(method == "W"), x = "n_of_seg", y = "count", facet.by = "cancer_type",
          xlab = FALSE, ylab = FALSE) +
  geom_abline(slope = 1, color = "red")
ggscatter(count_df2 %>% filter(method == "X"), x = "n_of_seg", y = "count", facet.by = "cancer_type",
          xlab = FALSE, ylab = FALSE) +
  geom_abline(slope = 1, color = "red")

count_summary <- count_df2 %>%
  group_by(cancer_type, sig_type, sample, method) %>%
  summarise(rmse = sqrt(sum((count - n_of_seg)^2) / length(count)),
            change_frac = 100 * (rmse / n_of_seg),
            .groups = "drop")

count_summary

p <- ggboxplot(count_summary %>%
            mutate(
              sig_type = as.integer(stringr::str_replace(sig_type, "_Sigs", ""))
            ),
          x = "sig_type", y = "rmse", fill = "method", size = 0.5, width = 0.5, outlier.size = 0.5,
          facet.by = "cancer_type", scales = "free_y",
          xlab = "Extracted signature number", ylab = "RMSE between observed segments and estimated segments") +
  stat_compare_means(aes(group = method, label = ..p.signif..))
ggsave("output/estimation_RMSE_comparison_between_methods.pdf", plot = p, width = 12, height = 7)

p <- ggboxplot(count_summary %>%
                 mutate(
                   sig_type = as.integer(stringr::str_replace(sig_type, "_Sigs", ""))
                 ),
               x = "sig_type", y = "change_frac", fill = "method", size = 0.5, width = 0.5, outlier.size = 0.5,
               facet.by = "cancer_type", scales = "free_y",
               xlab = "Extracted signature number", ylab = "Change between observed segments and estimated segments (%)") +
  stat_compare_means(aes(group = method, label = ..p.signif..))
ggsave("output/estimation_change_frac_comparison_between_methods.pdf", plot = p, width = 12, height = 7)
