library(sigminer)

pcawg_sig <- readRDS("data/pcawg_cn_sigs_CN176_signature.rds")

p <- show_sig_profile(pcawg_sig, style = "cosmic", mode = "copynumber", method = "X", font_scale = 0.7)
ggplot2::ggsave("output/pcawg_cn_sigs_c176.pdf", plot = p, width = 14, height = 8)

p <- show_sig_profile(pcawg_sig, style = "cosmic", mode = "copynumber", method = "X", font_scale = 0.7, by_context = TRUE, x_lab = NULL)
ggplot2::ggsave("output/pcawg_cn_sigs_c176_by_context.pdf", plot = p, width = 14, height = 8)

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
                           font_scale = 0.7, by_context = TRUE)

ggplot2::ggsave("output/pcawg_cn_sigs_loop_c176_by_context.pdf", plot = p, width = 14, height = 30, )
