# Cargar bibliotecas necesarias
rm(list = ls())
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


set.seed(1234567890, kind = "Mersenne-Twister")

sample.size <- 121
R <- 10000
mu <- 1
L <- 5
alpha1 <- -1.5
means <- numeric(R)
cvs <- numeric(R)

for (r in 1:R) {
  z <- gi0_sample(mu, alpha1, L, sample.size)
  means[r] <- mean(z)
  cvs[r] <- sd(z) / mean(z)
}


data <- data.frame(Mean = means, CV = cvs)

# 
start_values <- list(beta0 = 0.5, beta1 = 0.5)
model <- nls(CV ~ sqrt(sample.size) * (1 - exp(-(beta0 + beta1 * Mean))), 
             data = data, start = start_values)

# 
summary(model)

# 
ggplot(data, aes(x = Mean, y = CV)) +
  geom_point(alpha = 0.5, color = "red") +
  stat_function(fun = function(x) sqrt(sample.size) * (1 - exp(-(coef(model)[1] + coef(model)[2] * x))),
                color = "blue", size = 1) +
  labs(title = "Scatter plot of CV versus Mean for gI0 distribution",
       x = "Mean",
       y = "Coefficient of Variation (CV)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90")
  )