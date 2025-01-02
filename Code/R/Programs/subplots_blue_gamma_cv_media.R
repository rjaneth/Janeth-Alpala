# subplots_blue_gamma_cv_media
#rm(list = ls())
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

# Parámetros de la simulación
sample.size <- 49
R <- 10000
mu <- 1
L <- 5

# Valores específicos de alpha
alphas <- c(-1.5, -3, -5, -10, -50, -1000)

results <- data.frame()

for (alpha in alphas) {
  means <- numeric(R)
  cvs <- numeric(R)
  
  for (r in 1:R) {
    z <- gi0_sample(mu, alpha, L, sample.size)
    means[r] <- mean(z)
    cvs[r] <- sd(z) / mean(z)
  }
  
  temp_data <- data.frame(Mean = means, CV = cvs, Alpha = factor(alpha))
  results <- rbind(results, temp_data)
}

# Crear el gráfico con facet_wrap
p <- ggplot(results, aes(x = Mean, y = CV)) +
  geom_point(alpha = 0.5, color = "#56B4E9", size = 0.5) +
  labs(title = "",#Scatter Plot of CV versus Mean for Different Alpha Values
       x = "Mean",
       y = "CV") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90")
  ) +
  facet_wrap(~ Alpha, ncol = 3, scales = "free")

# Mostrar el gráfico
print(p)