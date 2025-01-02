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

# Regression Model Fitting
# Calculate y
data$y <- data$CV / sqrt(sample.size)

# Exclude data points where y >= 1 to avoid log of non-positive numbers
data_valid <- subset(data, y < 1)

# Calculate the response variable
data_valid$Response <- -log(1 - data_valid$y)

# Fit the linear model
model_lm <- lm(Response ~ Mean, data = data_valid)

# Display the summary of the model
summary(model_lm)

# Extract coefficients
beta0 <- coef(model_lm)[1]
beta1 <- coef(model_lm)[2]

# Compute fitted Response
data_valid$Fitted_Response <- predict(model_lm)

# Compute fitted y
data_valid$Fitted_y <- 1 - exp(-data_valid$Fitted_Response)

# Compute fitted CV
data_valid$Fitted_CV <- data_valid$Fitted_y * sqrt(sample.size)

# Plot the original data and fitted curve
plot(data_valid$Mean, data_valid$CV, xlab = "Mean", ylab = "CV",
     main = "CV vs Mean with Fitted Model", pch = 16, cex = 0.5)
lines(data_valid$Mean[order(data_valid$Mean)], 
      data_valid$Fitted_CV[order(data_valid$Mean)], 
      col = "blue", lwd = 2)
