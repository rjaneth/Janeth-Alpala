rm(list = ls())

library(ggplot2)
library(ggsci)
library(invgamma)
library(latex2exp)
library(data.table)

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

set.seed(1234567890, kind = "Mersenne-Twister")


alpha <- 0.05


R <- 100
mu <- 1
L <- 5
B <- 2
alpha1 <- -6


sample_sizes <- c(25, 49, 81, 121)


error_type_I_results <- numeric(length(sample_sizes))


for (i in seq_along(sample_sizes)) {
  s <- sample_sizes[i]
  false_positives_count <- 0
  
  for (r in 1:R) {
    #
    z0 <- gi0_sample(1, -6, 5, s)
    
    # prueba bajo H0
    TestStat <- sqrt(R) * (bootstrap_correa_estimator_log_mean(z0, B) + (-L + log(L) - lgamma(L) - (1 - L) * digamma(L))) / sd(z0)
    #TestStat <- bootstrap_correa_estimator_log_mean(z0, B) + (-L + log(L) - lgamma(L) - (1 - L) * digamma(L)) 
    

    
    # 
    if (abs(TestStat) > qnorm(1 - alpha/2)) {
      false_positives_count <- false_positives_count + 1
    }
  }
  
 
  cat(paste("False Positives for Sample Size", s, ":", false_positives_count, "\n"))
  
  # Calcular y almacenar el error tipo I
  error_type_I <- false_positives_count / R
  error_type_I_results[i] <- error_type_I
  
  cat(paste("Error Type I for Sample Size", s, ":", error_type_I, "\n"))
}

# Graficar resultados
plot(sample_sizes, error_type_I_results, type = "b", pch = 16, col = "red",
     xlab = "Sample Size", ylab = "Error Type I",
     main = "Error Type I Analysis")
