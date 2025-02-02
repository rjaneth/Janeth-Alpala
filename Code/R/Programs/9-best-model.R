# Load necessary libraries
# Install 'minpack.lm' package if you haven't already
# install.packages("minpack.lm")
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
R <- 10000
mu <- 1
L <- 5
alpha1 <- -1.5

# Initialize vectors to store means and CVs
means <- numeric(R)
cvs <- numeric(R)

# Generate samples and compute means and CVs
for (r in 1:R) {
  z <- gi0_sample(mu, alpha1, L, sample.size)
  means[r] <- mean(z)
  cvs[r] <- sd(z) / mean(z)
}

# Create data frame
data <- data.frame(Mean = means, CV = cvs)

library(minpack.lm)

# Use nlsLM for robust fitting
nls_model_exp_simplified_lm <- nlsLM(CV ~ sqrt(sample.size) * (1 - exp(-(beta0 + beta1 * Mean)) - beta2),
                                     data = data, start = list(beta0 = 0.1, beta1 = 0.01, beta2 = 0.01))

summary(nls_model_exp_simplified_lm)

data$Fitted_CV_simplified_lm <- predict(nls_model_exp_simplified_lm)

plot(data$Mean, data$CV, xlab = "Mean", ylab = "CV", main = paste("CV vs Mean, alpha =", alpha1), pch = 16, cex = 0.5)
lines(data$Mean[order(data$Mean)], data$Fitted_CV_simplified_lm[order(data$Mean)], col = "blue", lwd = 2)
