

rm(list = ls())

library(ggplot2)
library(ggsci)
library(invgamma)
library(latex2exp)
library(data.table)
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




alpha <- 0.05


R <- 3000
mu <- 1
L <- 5
B <- 2
alpha1 <- -600

sample_sizes <- c( 25, 49, 81, 121)


power_results <- numeric(length(sample_sizes))

# Escenarios bajo H1 (por ejemplo, clutter heterogéneo)
for (i in seq_along(sample_sizes)) {
  s <- sample_sizes[i]
  reject_H0_count <- 0
  
  for (r in 1:R) {
   
    z1 <- gi0_sample(1, -600, 5, s)
    
   
    TestStat <- sqrt(R) * (bootstrap_correa_estimator_log_mean(z1, B) + (-L + log(L) - lgamma(L) - (1 - L) * digamma(L))) / sd(z1)
    
    
    if (abs(TestStat) > qnorm(1 - alpha/2)) {
      reject_H0_count <- reject_H0_count + 1
    }
  }
  
  # Calcular y almacenar el poder (1 - beta)
  power <- reject_H0_count / R
  power_results[i] <- power
  
  cat(paste("Power (1 - beta) for Sample Size", s, ":", power, "\n"))
}

# Graficar resultados
plot(sample_sizes, power_results, type = "b", pch = 16, col = "blue",
     xlab = "Sample Size", ylab = "Power (1 - beta)",
     main = "Power Analysis")
