# Limpiar el entorno
rm(list = ls())
library(ggplot2)
library(ggsci)
library(univariateML)
library(compiler)
library(MASS)
library(fitdistrplus)
library(invgamma)
library(gridExtra)
# Cargar las librerías necesarias
library(ggplot2)
library(gganimate)
library(ggsci)
library(magick)

source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/entropy_gamma_sar.R")
source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")


set.seed(1234567890, kind = "Mersenne-Twister")


sample.size <- 49
R <- 1000
mu <- 1
L <- 5
alphas <- seq(-1.5, -20, by = -1.0)


results <- data.frame()


for (alpha in alphas) {
  means <- numeric(R)
  cvs <- numeric(R)
  
  for (r in 1:R) {
    z <- gi0_sample(mu, alpha, L, sample.size)
    means[r] <- mean(z)
    cvs[r] <- sd(z) / mean(z)
  }
  
  
  temp_data <- data.frame(Mean = means, CV = cvs, Alpha = alpha)
  results <- rbind(results, temp_data)
}

# 
p <- ggplot(results, aes(x = Mean, y = CV)) +
  geom_point(alpha = 0.5, color = "red") +
  labs(title = "Scatter plot of CV versus Mean for gI0 distribution",
       x = "Mean",
       y = "Coefficient of Variation (CV)",
       subtitle = "Alpha: {closest_state}") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90")
  ) +
  transition_states(Alpha, transition_length = 2, state_length = 1) +
  ease_aes('linear')

#  
animate(p, nframes = 100, width = 800, height = 600, renderer = gifski_renderer())
