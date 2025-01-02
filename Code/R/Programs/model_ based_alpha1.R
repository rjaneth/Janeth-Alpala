#533 modelo basado en alpha1
rm(list = ls())
library(minpack.lm)
library(invgamma)
source("../MainFunctions/gi0_sample.R")

# Asegúrate de que la función gi0_sample está disponible
# source("../MainFunctions/gi0_sample.R")

# Definimos la función fit_cv_mean_model
fit_cv_mean_model <- function(mu, alpha1, L, sample.size, R) {
  # Generar datos
  means <- numeric(R)
  cvs <- numeric(R)
  
  for (r in 1:R) {
    z <- gi0_sample(mu, alpha1, L, sample.size)
    means[r] <- mean(z)
    cvs[r] <- sd(z) / mean(z)
  }
  
  data <- data.frame(Mean = means, CV = cvs)
  
  # Determinar el modelo basado en alpha1
  if (alpha1 >= -1.8) {
    # Usar el modelo no lineal
    nls_model <- nlsLM(CV ~ sqrt(sample.size) * (1 - exp(-(beta0 + beta1 * Mean)) - beta2),
                       data = data, start = list(beta0 = 0.01, beta1 = 0.01, beta2 = 0.01),
                       control = nls.lm.control(maxiter = 500))
    model_type <- "No lineal"
  } else {
    # Usar un modelo lineal
    nls_model <- lm(CV ~ Mean, data = data)
    model_type <- "Lineal"
  }
  
  # Salida
  list(model = nls_model, data = data, model_type = model_type)
}
#set.seed(1234567890, kind = "Mersenne-Twister")
# Parámetros
mu <- 1
alpha1 <- -3.5
L <- 5
sample.size <- 49
R <- 30000

# Ajustar el modelo
resultado <- fit_cv_mean_model(mu, alpha1, L, sample.size, R)

# Mostrar resumen del modelo
summary(resultado$model)

# Graficar los datos
plot(resultado$data$Mean, resultado$data$CV, xlab = "Mean", ylab = "CV",
     main = paste("CV vs Mean, alpha =", alpha1), pch = 16, cex = 0.7)

# Agregar la línea de ajuste del modelo
if (resultado$model_type == "No lineal") {
  # Obtener predicciones del modelo no lineal
  resultado$data$Fitted_CV <- predict(resultado$model)
  # Ordenar los datos para una línea suave
  ord <- order(resultado$data$Mean)
  lines(resultado$data$Mean[ord], resultado$data$Fitted_CV[ord], col = "blue", lwd = 2)
} else {
  # Modelo lineal
  abline(resultado$model, col = "red", lwd = 2)
}

# Leyenda para indicar el tipo de modelo
legend("topright", legend = paste("Modelo", resultado$model_type), col = ifelse(resultado$model_type == "No lineal", "blue", "red"), lwd = 2)
