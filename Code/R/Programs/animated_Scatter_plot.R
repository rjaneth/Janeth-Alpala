# Limpiar el entorno
rm(list = ls())

# Cargar las librerías necesarias
library(ggplot2)
library(ggsci)
library(gganimate)
library(invgamma)

source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/entropy_gamma_sar.R")
source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")

# Fijar la semilla para reproducibilidad
set.seed(1234567890, kind = "Mersenne-Twister")

# Definir parámetros
sample.size <- 49
R <- 3000
mu <- 1
L <- 5
alphas <- seq(-1.5, -10, by = -0.5)

# Inicializar data frame para almacenar los resultados
results <- data.frame()

# Generar las muestras y calcular la media y el CV para cada alpha
for (alpha in alphas) {
  means <- numeric(R)
  cvs <- numeric(R)
  
  for (r in 1:R) {
    z <- gi0_sample(mu, alpha, L, sample.size)
    means[r] <- mean(z)
    cvs[r] <- sd(z) / mean(z)
  }
  
  # Agregar los resultados al data frame
  temp_data <- data.frame(Mean = means, CV = cvs, Alpha = alpha)
  results <- rbind(results, temp_data)
}

# Crear el plot animado
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
  transition_states(Alpha, transition_length = 2, state_length = 1)

# Guardar la animación
#animate(p, nframes = 100, width = 800, height = 600)
#anim_save("gi0_animation.gif")
