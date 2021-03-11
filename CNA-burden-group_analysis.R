# 先按 cancer_type 画分布即可

library(ggridges)

tidy_info <- readRDS("data/pcawg_sample_tidy_info.rds")

p <- ggplot(tidy_info, aes(x = cna_burden, y = cancer_type, fill = cancer_type)) +
  geom_density_ridges() +
  labs(x = "CNA burden", y = NULL) +
  theme_ridges() +
  theme(legend.position = "none")

ggsave("output/cna_burden_pancan_dist.pdf", plot = p, width = 6, height = 7)

# Split samples into groups based on CNA burden distribution

# library(flexmix)
#
# N1 <- 100
# N2 <- 10
#
# a <- rpois(N1, 0)
# b <- rpois(N2, 50)
#
# x <- c(a,b)
# class <- c(rep('a', N1), rep('b', N2))
# data <- data.frame(cbind(x=as.numeric(x), class=as.factor(class)))
#
# fit3 <- flexmix(x ~ 1, data = data, k = 2, model = flexmix::FLXMCmvpois())
# fit3 <- flexmix(x ~ 1, data = data, k = 2, model = FLXMCmvnorm())
# clusters(fit3)
