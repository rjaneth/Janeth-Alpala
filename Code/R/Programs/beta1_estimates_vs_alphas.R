# 530 beta1_estimates vs alphas
# Initialize vectors to store parameter estimates

library(minpack.lm)

# Cargar bibliotecas necesarias
#rm(list = ls())
library(ggplot2)
#library(nls2)
library(ggsci)
library(univariateML)
library(compiler)
library(MASS)
library(fitdistrplus)
library(invgamma)
library(gridExtra)

library(ggplot2)
library(ggsci)


source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/entropy_gamma_sar.R")
source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")

# Set parameters
set.seed(1234567890, kind = "Mersenne-Twister")
sample.size <- 49

mu <- 1

alphas <- seq(-1.2, -2.8, by = -0.1)
beta0_estimates <- numeric(length(alphas))
beta1_estimates <- numeric(length(alphas))
beta2_estimates <- numeric(length(alphas))

for (i in seq_along(alphas)) {
  alpha1 <- alphas[i]
  
  # Generate data (use a smaller R for quicker execution)
  R <- 300
  means <- numeric(R)
  cvs <- numeric(R)
  
  for (r in 1:R) {
    z <- gi0_sample(mu, alpha1, L, sample.size)
    means[r] <- mean(z)
    cvs[r] <- sd(z) / mean(z)
  }
  
  data <- data.frame(Mean = means, CV = cvs)
  
  # Fit the model with initial estimates based on previous iteration
  if (i == 1) {
    # Initial starting values
    start_list <- list(beta0 = 0.01, beta1 = 0.01, beta2 = 0.01)
  } else {
    # Use previous estimates as starting values
    start_list <- list(beta0 = beta0_estimates[i - 1], 
                       beta1 = beta1_estimates[i - 1], 
                       beta2 = beta2_estimates[i - 1])
  }
  
  nls_model <- try(nlsLM(CV ~ sqrt(sample.size) * (1 - exp(-(beta0 + beta1 * Mean)) - beta2),
                         data = data, start = start_list, 
                         control = nls.lm.control(maxiter = 500)), silent = TRUE)
  
  if (class(nls_model) != "try-error") {
    coef_estimates <- coef(nls_model)
    beta0_estimates[i] <- coef_estimates["beta0"]
    beta1_estimates[i] <- coef_estimates["beta1"]
    beta2_estimates[i] <- coef_estimates["beta2"]
  } else {
    beta0_estimates[i] <- NA
    beta1_estimates[i] <- NA
    beta2_estimates[i] <- NA
  }
}

# Plot the parameter estimates against alpha
plot(alphas, beta1_estimates, type = "b", xlab = "Alpha", ylab = "Beta1 Estimate",
     main = "Beta1 Estimates vs Alpha")

# # Fit a linear model to beta1_estimates vs alphas (excluding NA values)
# valid_indices <- which(!is.na(beta1_estimates))
# lm_beta1 <- lm(beta1_estimates[valid_indices] ~ alphas[valid_indices])
# 
# # Predict starting beta1 for new alpha
# alpha_new <- -3.0
# beta1_start <- predict(lm_beta1, newdata = data.frame(alphas = alpha_new))
# 
# data$LogMean <- log(data$Mean)
# data$LogCV <- log(data$CV)
# 
# # Fit a linear model
# lm_model <- lm(LogCV ~ LogMean, data = data)
# summary(lm_model)
