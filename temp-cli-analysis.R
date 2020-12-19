library(sigminer)
library(tidyverse)

pcawg_types <- readRDS("data/pcawg_type_info.rds")
pcawg_activity <- readRDS("data/pcawg_cn_sigs_CN176_activity.rds")
cli <- readRDS("data/pcawg_samp_info_sp.rds")
os <- cli$pcawg_donor_clinical_August2016_v9_sp %>%
  select(xena_sample, donor_vital_status, donor_survival_time) %>%
  set_names(c("sample", "os", "time")) %>%
  na.omit() %>%
  mutate(os = ifelse(os == "deceased", 1, 0))

keep_samps <- pcawg_activity$similarity >= 0.75
df <- merge(pcawg_activity$abs_activity[keep_samps], pcawg_types, by = "sample")
df <- left_join(df, os, by = "sample")

df2 <- df %>%
  mutate_at(vars(starts_with("Sig")), ~ . / 10)

library(ezcox)

p1 <- show_forest(df,
  covariates = paste0("Sig", 1:11),
  time = "time", status = "os",
  merge_models = TRUE, add_caption = FALSE, point_size = 2
)

p2 <- show_forest(df,
  covariates = "Sig1", controls = paste0("Sig", 2:11),
  time = "time", status = "os",
  merge_models = TRUE, add_caption = FALSE, point_size = 2
)

dir.create("output/pcawg_os_cox")
ggsave("output/pcawg_os_cox/raw_activity_unicox.pdf", plot = p1, width = 7, height = 5)
ggsave("output/pcawg_os_cox/raw_activity_multicox.pdf", plot = p2, width = 7, height = 5)

p3 <- show_forest(df2,
  covariates = paste0("Sig", 1:11),
  time = "time", status = "os",
  merge_models = TRUE, add_caption = FALSE, point_size = 2
)

p4 <- show_forest(df2,
  covariates = "Sig1", controls = paste0("Sig", 2:11),
  time = "time", status = "os",
  merge_models = TRUE, add_caption = FALSE, point_size = 2
)

ggsave("output/pcawg_os_cox/raw_activity_div10_unicox.pdf", plot = p3, width = 7, height = 5)
ggsave("output/pcawg_os_cox/raw_activity_div10_multicox.pdf", plot = p4, width = 7, height = 5)

df3 <- df %>%
  mutate_at(vars(starts_with("Sig")), ~ ifelse(. > median(., na.rm = TRUE), 1, 0))
df3

summary(df3)

p5 <- show_forest(df3,
  covariates = paste0("Sig", 1:11),
  time = "time", status = "os",
  merge_models = TRUE, add_caption = FALSE, point_size = 2
)

p6 <- show_forest(df3,
  covariates = "Sig1", controls = paste0("Sig", 2:11),
  time = "time", status = "os",
  merge_models = TRUE, add_caption = FALSE, point_size = 2
)

ggsave("output/pcawg_os_cox/two_grp_activity_unicox.pdf", plot = p5, width = 7, height = 5)
ggsave("output/pcawg_os_cox/two_grp_activity_multicox.pdf", plot = p6, width = 7, height = 5)


# Redo for cancer types ---------------------------------------------------

table(na.omit(df)$cancer_type)

typeList <- list(
  Eso_AdenoCA = c("Eso-AdenoCA"),
  Lung = c("Lung-AdenoCA", "Lung-SCC"),
  Ovary_AdenoCA = c("Ovary-AdenoCA"),
  Panc_AdenoCA = c("Panc-AdenoCA"),
  Prost_AdenoCA = c("Prost-AdenoCA"),
  SKCM = c("Skin-Melanoma")
)

for (i in names(typeList)) {
  type = typeList[[i]]

  p3 <- show_forest(df2 %>% filter(cancer_type %in% type),
                    covariates = paste0("Sig", 1:11),
                    time = "time", status = "os",
                    merge_models = TRUE, add_caption = FALSE, point_size = 2
  )

  p4 <- show_forest(df2 %>% filter(cancer_type %in% type),
                    covariates = "Sig1", controls = paste0("Sig", 2:11),
                    time = "time", status = "os",
                    merge_models = TRUE, add_caption = FALSE, point_size = 2
  )

  ggsave(paste0("output/pcawg_os_cox/raw_activity_div10_unicox_", i, ".pdf"),
         plot = p3, width = 7, height = 5)
  ggsave(paste0("output/pcawg_os_cox/raw_activity_div10_multicox_", i, ".pdf"),
         plot = p4, width = 7, height = 5)

  p5 <- show_forest(df3 %>% filter(cancer_type %in% type),
                    covariates = paste0("Sig", 1:11),
                    time = "time", status = "os",
                    merge_models = TRUE, add_caption = FALSE, point_size = 2
  )

  p6 <- show_forest(df3 %>% filter(cancer_type %in% type),
                    covariates = "Sig1", controls = paste0("Sig", 2:11),
                    time = "time", status = "os",
                    merge_models = TRUE, add_caption = FALSE, point_size = 2
  )

  ggsave(paste0("output/pcawg_os_cox/two_grp_activity_unicox_", i, ".pdf"),
         plot = p5, width = 7, height = 5)
  ggsave(paste0("output/pcawg_os_cox/two_grp_activity_multicox_", i, ".pdf"),
         plot = p6, width = 7, height = 5)
}

