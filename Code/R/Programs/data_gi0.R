rm(list = ls())

library(ggplot2)
library(ggsci)
library(invgamma)
library(latex2exp)
if(!require("rstudioapi")) install("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../MainFunctions/ebrahimi_estimator.R")
source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/correa_estimator.R")
source("../MainFunctions/bootstrap_correa_estimator.R")
source("../MainFunctions/bootstrap_correa_estimator_log_mean.R")

source("../MainFunctions/al_omari_1_estimator.R")
source("../MainFunctions/bootstrap_al_omari_1_estimator.R")

source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")

set.seed(1234567890, kind="Mersenne-Twister")

R <- 1000
mu <- 1
L <- 5
B <- 2
alpha1<- -600
sample.size <- c(9, 25, 49, 81)



TestStatistics <- NULL

for(s in sample.size){
  for(r in 1:R){
    z <- gi0_sample(mu, alpha1, L,  s)
   
    #TestStat <- log(mean(z))  - (-L + log(L) - lgamma(L) - (1 - L) * digamma(L))
    TestStat <- bootstrap_correa_estimator_log_mean(z, B) + (-L + log(L) - lgamma(L) - (1 - L) * digamma(L))
    #TestStat <- bootstrap_correa_estimator(z, B) -log(mean(z))  + (-L + log(L) - lgamma(L) - (1 - L) * digamma(L)+L + lgamma(L-alpha1) - (L-alpha1)*(digamma(L - alpha1))+(1-alpha1)*digamma(- alpha1)-log(-1 - alpha1)-lgamma(-alpha1) ) 
    TestStatistics <- rbind(TestStatistics,
                            c(s, TestStat))
  }
  mean_entropy1 <- mean(TestStat)
  cat(" mean", mean_entropy1, "\n" )
}

TestStatistics <- data.frame(TestStatistics)
names(TestStatistics) <- c("Sample Size", "Test Statistics")




ggplot(TestStatistics, aes(x = `Test Statistics`, col = factor(`Sample Size`), linetype = factor(`Sample Size`))) +
  geom_line(stat = "density", linewidth = 1.0) +  
  scale_color_manual(
    values = pal_cosmic()(length(sample.size)),
    name = "Sample Size"
  ) +
  scale_linetype_manual(
    values = rep("solid", length(sample.size)),  
    name = "Sample Size"
  ) +
  
  labs(x = "Test Statistics", y = "Density") +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    legend.position = "bottom",
    legend.key.size = unit(1, "lines")  
  )
