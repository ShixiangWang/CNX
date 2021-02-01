library(sigminer)

pcawg_sig <- readRDS("data/pcawg_cn_sigs_CN176_signature.rds")

p <- show_sig_profile(pcawg_sig, style = "cosmic", mode = "copynumber", method = "X", font_scale = 0.7)
ggplot2::ggsave("output/pcawg_cn_sigs_c176.pdf", plot = p, width = 14, height = 8)

p <- show_sig_profile(pcawg_sig, style = "cosmic", mode = "copynumber", method = "X", font_scale = 0.7, by_context = TRUE, x_lab = NULL)
ggplot2::ggsave("output/pcawg_cn_sigs_c176_by_context.pdf", plot = p, width = 14, height = 10)

p2 <- show_sig_profile(pcawg_sig, style = "cosmic", mode = "copynumber", method = "X", font_scale = 0.6,
                       rm_axis_text = TRUE, x_lab = NULL)
ggplot2::ggsave("output/pcawg_cn_sigs_c176_without_axis_text.pdf", plot = p2, width = 12, height = 7)


p2 <- show_sig_profile(pcawg_sig, style = "cosmic", mode = "copynumber", method = "X", font_scale = 0.6,
                       rm_axis_text = TRUE, x_lab = NULL, by_context = TRUE)
ggplot2::ggsave("output/pcawg_cn_sigs_c176_without_axis_text_by_context.pdf", plot = p2, width = 12, height = 7)

p <- show_sig_profile_loop(pcawg_sig,
                           style = "cosmic",
                           mode = "copynumber",
                           method = "X",
                           font_scale = 0.7)

ggplot2::ggsave("output/pcawg_cn_sigs_loop_c176.pdf", plot = p, width = 14, height = 22)


p <- show_sig_profile_loop(pcawg_sig,
                           x_lab = NULL,
                           style = "cosmic",
                           mode = "copynumber",
                           method = "X",
                           font_scale = 0.6, by_context = TRUE)

ggplot2::ggsave("output/pcawg_cn_sigs_loop_c176_by_context.pdf", plot = p, width = 14, height = 30)


tcga_sigs <- readRDS("data/TCGA/tcga_cn_sigs_CN176_signature.rds")
sapply(get_sig_exposure(tcga_sigs)[, -1], sum)

p <- show_sig_profile_loop(tcga_sigs,
                           x_lab = NULL,
                           style = "cosmic",
                           mode = "copynumber",
                           method = "X",
                           font_scale = 0.6, by_context = TRUE)

ggplot2::ggsave("output/tcga_cn_sigs_loop_c176_by_context.pdf", plot = p, width = 14, height = 40)


# One representative sample -----------------------------------------------

s <- "SP112845"
obj <- readRDS("data/pcawg_cn_obj.rds")
tally_x <- readRDS("data/pcawg_cn_tally_X.rds")

show_cn_profile(obj, sample = s)
p <- show_catalogue(t(tally_x$nmf_matrix), samples = s,
               style = "cosmic", mode = "copynumber", method = "X",
               by_context = TRUE, font_scale = 0.7, normalize = "row")
ggplot2::ggsave("output/showcase_catalog_profile.pdf", plot = p, width = 14, height = 2.5)

reconstructed_mat <- pcawg_sig$Signature.norm %*% pcawg_sig$Exposure
p <- show_catalogue(reconstructed_mat, samples = s,
                    style = "cosmic", mode = "copynumber", method = "X",
                    by_context = TRUE, font_scale = 0.7, normalize = "row")
ggplot2::ggsave("output/showcase_rect_catalog_profile.pdf", plot = p, width = 14, height = 2.5)

cosine(reconstructed_mat[, s], t(tally_x$nmf_matrix)[, s])

act <- readRDS("data/pcawg_cn_sigs_CN176_activity.rds")
round(act$relative[sample == s][, -1], 3)

p <- show_sig_profile(pcawg_sig, sig_names = paste0("CNS", c(3, 8, 10, 11, 14)),
                      style = "cosmic", mode = "copynumber", method = "X",
                      by_context = TRUE,  font_scale = 0.7)
ggplot2::ggsave("output/showcase_sig_profile.pdf", plot = p, width = 14, height = 5)

x = get_sig_exposure(pcawg_sig, type = "relative")
round(x[sample == "SP112845"][,-1], 3)
