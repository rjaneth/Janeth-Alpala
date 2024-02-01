
rm(list = ls())

library(ggplot2)
library(ggsci)
library(invgamma)
library(latex2exp)
library(data.table)

# if(!require("rstudioapi")) install("rstudioapi")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../MainFunctions/ebrahimi_estimator.R")
source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/correa_estimator.R")
source("../MainFunctions/bootstrap_correa_estimator.R")
#source("../MainFunctions/bootstrap_correa_estimator_log_mean.R")

source("../MainFunctions/al_omari_1_estimator.R")
source("../MainFunctions/bootstrap_al_omari_1_estimator.R")

source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")
source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/bootstrap_correa_estimator_log_mean.R")
set.seed(1234567890, kind = "Mersenne-Twister")

# Function to calculate p-value
calculate_p_value <- function(test_statistic, mu_W, sigma_W) {
  epsilon <- test_statistic  / sigma_W
  p_value <- 2 * (1 - pnorm(abs(epsilon)))
  
  #cat("epsilon:", epsilon, "\n")
  return(p_value)
}

R <- 1000
mu <- 1
L <- 5
B <- 1
sample.size <- c( 50, 100)



TestStatistics1 <- NULL

mean_entropy1 <- numeric(length(sample.size))
sd_entropy1 <- numeric(length(sample.size))

for (s in sample.size) {
  TestStat1 <- numeric(R)
  
  for (r in 1:R) {
    z <- gamma_sar_sample(L, mu, s) #H0: los datos provienen de gamma sar
    #z <- gi0_sample(mu, alpha1, L, s) #H1: los datos provienen de gi0
    TestStat1[r] <- bootstrap_correa_estimator_log_mean(z, B) - (L - log(L) + lgamma(L) + (1 - L) * digamma(L)
    )
  }
  mean_entropy1[sample.size == s] <- mean(TestStat1)
  sd_entropy1[sample.size == s] <- sd(TestStat1)
  
  TestStatistics1 <- rbind(TestStatistics1, data.frame("Sample_Size" = rep(s, R), "Test_Statistics" = TestStat1))
}

# Calculate p-values
mu_W <-  mean_entropy1  
sigma_W <- sqrt(sd_entropy1^2)  

# Calculate p-values
p_values <- apply(TestStatistics1, 1, function(row) {
  calculate_p_value(row["Test_Statistics"], mu_W, sigma_W[row["Sample_Size"] == s])
})


# 
result <- data.frame("Sample_Size" = TestStatistics1$Sample_Size, "P_Value" = p_values)

#print(result)


# 
alpha_nominal <- 0.05

#  tasa de error de Tipo I
type_I_error_rate <- sum(result$P_Value < alpha_nominal) / nrow(result)

# 
cat("Type I Error Rate:", type_I_error_rate, "\n")
