# Limpiar el entorno
#rm(list = ls())

# Cargar librerías necesarias
library(ggplot2)
library(ggsci)
library(invgamma)

# Cargar funciones necesarias
source("../MainFunctions/ebrahimi_estimator.R")
source("../MainFunctions/gamma_sar_sample.R") 
source("../MainFunctions/renyi_entropy_estimator_v1.R")
source("../MainFunctions/renyi_entropy_estimator.R")
source("../MainFunctions/correa_estimator.R")
source("../MainFunctions/bootstrap_renyi_entropy_estimator_v1.R")
source("../MainFunctions/bootstrap_renyi_entropy_estimator.R")
source("../MainFunctions/bootstrap_ebrahimi_estimator.R")
source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")

# Configuración de la semilla
set.seed(1234567890, kind = "Mersenne-Twister")

# Parámetros de simulación
R <- 200
mu <- 100
L <- 5
B <- 10
lambda <- 0.8
sample.size <-81

# Para medir el tiempo de ejecución de todo el proceso
total_time <- system.time({
  # Para guardar los resultados
  TestStatistics <- numeric(R)
  
  for (r in 1:R) {
    # Generar muestra gamma SAR
    z <- gamma_sar_sample(L, mu, sample.size)
    
    # Calcular el estadístico usando el bootstrap del estimador de Rényi
    TestStat <- bootstrap_renyi_entropy_estimator_v1(z, B, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
                                                                        (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(z))-log(L))
    
    
    TestStatistics[r] <- TestStat
  }
  
  # Calcular la entropía media
  mean_entropy <- mean(TestStatistics, na.rm = TRUE)
  cat("Mean Entropy:", mean_entropy, "\n")
})

# Imprimir el tiempo total de ejecución
print(total_time)

# Convertir resultados en un data frame para la visualización
TestStatistics_df <- data.frame(
  Sample_Size = rep(sample.size, R),
  Test_Statistics = TestStatistics
)

# Visualización de los resultados
ggplot(TestStatistics_df, aes(x = Test_Statistics, col = factor(Sample_Size), linetype = factor(Sample_Size))) +
  geom_line(stat = "density", linewidth = 1.0) +  
  scale_color_manual(
    values = pal_cosmic()(length(unique(TestStatistics_df$Sample_Size))),
    name = "Sample Size"
  ) +
  scale_linetype_manual(
    values = rep("solid", length(unique(TestStatistics_df$Sample_Size))),
    name = "Sample Size"
  ) +
  labs(x = "Test Statistics", y = "Density") +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    legend.position = "bottom",
    legend.key.size = unit(1, "lines")
  )