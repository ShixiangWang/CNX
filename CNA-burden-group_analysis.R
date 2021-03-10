# Split samples into groups based on CNA burden distribution

library(flexmix)

N1 <- 100
N2 <- 10

a <- rpois(N1, 0)
b <- rpois(N2, 50)

x <- c(a,b)
class <- c(rep('a', N1), rep('b', N2))
data <- data.frame(cbind(x=as.numeric(x), class=as.factor(class)))

fit3 <- flexmix(x ~ 1, data = data, k = 2, model = flexmix::FLXMCmvpois())
fit3 <- flexmix(x ~ 1, data = data, k = 2, model = FLXMCmvnorm())
clusters(fit3)
