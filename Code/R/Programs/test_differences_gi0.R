rm(list = ls())

library(ggplot2)
library(ggsci)
library(invgamma)
library(latex2exp)
source("../MainFunctions/ebrahimi_estimator.R")
source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/correa_estimator.R")
source("../MainFunctions/bootstrap_correa_estimator.R")
source("../MainFunctions/bootstrap_correa_estimator_log_mean.R")

source("../MainFunctions/al_omari_1_estimator.R")
source("../MainFunctions/bootstrap_al_omari_1_estimator.R")

source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")

set.seed(1234567890, kind = "Mersenne-Twister")

R <- 1000
mu <- 1
L <- 5
B <- 2
alpha1 <- -6
sample.size <- c(81)

TestStatistics1 <- NULL
TestStatistics2 <- NULL
TestStatistics3 <- NULL
TestStatistics4 <- NULL
mean_entropy1 <- numeric(length(sample.size))
mean_entropy2 <- numeric(length(sample.size))
mean_entropy3 <- numeric(length(sample.size))
mean_entropy4 <- numeric(length(sample.size))
for (s in sample.size) {
  TestStat1 <- numeric(R)
  TestStat2 <- numeric(R)
  TestStat3 <- numeric(R)
  TestStat4 <- numeric(R)

  
  for (r in 1:R) {
    z <- gi0_sample(mu, alpha1, L, s)
    TestStat1[r] <- mean(z)  #- (-L + log(L) - lgamma(L) - (1 - L) * digamma(L))
    TestStat2[r] <- bootstrap_correa_estimator_log_mean(z, B) + (-L + log(L) - lgamma(L) - (1 - L) * digamma(L) + L + lgamma(L - alpha1) - (L - alpha1) * (digamma(L - alpha1)) + (1 - alpha1) * digamma(-alpha1) - log(-1 - alpha1) - lgamma(-alpha1))
    TestStat3[r] <- bootstrap_correa_estimator_log_mean(z, B) + (-L + log(L) - lgamma(L) - (1 - L) * digamma(L))
    TestStat4[r] <- bootstrap_correa_estimator(z, B)
  }
  
  mean_entropy1[sample.size == s] <- mean(TestStat1)
  cat(" mean for TestStat1", mean_entropy1[sample.size == s], "\n")
  
  mean_entropy2[sample.size == s] <- mean(TestStat2)
  cat(" mean for TestStat2", mean_entropy2[sample.size == s], "\n")
  mean_entropy3[sample.size == s] <- mean(TestStat3)
  cat(" mean for TestStat3", mean_entropy3[sample.size == s], "\n")
  
  mean_entropy4[sample.size == s] <- mean(TestStat4)
  cat(" mean for TestStat4", mean_entropy4[sample.size == s], "\n")
  
  
  TestStatistics1 <- rbind(TestStatistics1, data.frame("Sample_Size" = rep(s, R), "Test_Statistics" = TestStat1))
  TestStatistics2 <- rbind(TestStatistics2, data.frame("Sample_Size2" = rep(s, R), "Test_Statistics2" = TestStat2))
  TestStatistics3 <- rbind(TestStatistics3, data.frame("Sample_Size3" = rep(s, R), "Test_Statistics3" = TestStat3))
  TestStatistics4 <- rbind(TestStatistics4, data.frame("Sample_Size4" = rep(s, R), "Test_Statistics4" = TestStat4))
}

difference.betweeen.GammaSAR.GI0.n.Look <- function(alpha, L) {
  return(
    -L - lgamma(L-alpha) + 
      (L-alpha)*(digamma(L - alpha))- 
      (1-alpha)*digamma(- alpha)+
      log(-1 - alpha)+
      lgamma(-alpha)
  )
}

max_density1 <- max(density(TestStatistics1$Test_Statistics)$y, na.rm = TRUE)
max_density2 <- max(density(TestStatistics2$Test_Statistics2)$y, na.rm = TRUE)
#max_density3 <- max(density(TestStatistics3$Test_Statistics3)$y, na.rm = TRUE)

ggplot(TestStatistics1, aes(x = Test_Statistics, col = factor(Sample_Size), linetype = factor(Sample_Size))) +
  geom_line(stat = "density", linewidth = 1.0) +
  scale_color_manual(values = rainbow(7)[1], name = "Sample Size") +
  scale_linetype_manual(values = rep("solid", length(sample.size)), name = "Sample Size") +
  geom_line(data = TestStatistics2, aes(x = Test_Statistics2, col = factor(Sample_Size2), linetype = factor(Sample_Size2)), stat = "density", color = "blue", linetype = "dashed", linewidth = 1) +
  geom_line(data = TestStatistics3, aes(x = Test_Statistics3, col = factor(Sample_Size3), linetype = factor(Sample_Size3)), stat = "density", color = "green", linetype = "dashed", linewidth = 1) +
  geom_line(data = TestStatistics4, aes(x = Test_Statistics4, col = factor(Sample_Size4), linetype = factor(Sample_Size4)), stat = "density", color = "black", linetype = "dashed", linewidth = 1) +
  annotate("text", x = max(TestStatistics1$Test_Statistics) + 0.1, y = max_density1,
           label = TeX("$\\Gamma_{\\tiny{SAR}}$"), vjust = 5, hjust = 5, color = "red") +
  annotate("text", x = max(TestStatistics2$Test_Statistics2) + 0.1, y = max_density2,
           label = TeX("$G_I^0$"), vjust = 5, hjust = 5, color = "blue") +
  geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.n.Look(-600, 5), col="red") +
  annotate("text", x=difference.betweeen.GammaSAR.GI0.n.Look(-600, 5), y=10, parse=TRUE, label="alpha==-6", hjust=-.1, vjust=-1, col="red") +
  labs(x = "Test Statistics", y = "Density") +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    legend.position = "bottom",
    legend.key.size = unit(1, "lines")
  )
