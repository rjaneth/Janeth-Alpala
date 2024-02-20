
# rm(list = ls())
# if(!require("rstudioapi")) install("rstudioapi")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(e1071)
library(nortest)
library(fitdistrplus)

#load("./Data/results_data_Flevoland_cv_300_5.Rdata")
load("./Data/results_data_Flevoland_300_9.Rdata")
#summary(difference_values)


#load("results_data_Flevoland_300_9.Rdata")


p_values_matrix <- matrix(NA, nrow = nrow(difference_values), ncol = ncol(difference_values))


for (i in 1:nrow(difference_values)) {
  for (j in 1:ncol(difference_values)) {
    # Calcular z-score para cada diferencia
    z_score <- (difference_values[i, j] - mean(difference_values)) / sd(difference_values)
    
    # Calcular p-value para una prueba de dos colas
    p_value <- 2 * pnorm(-abs(z_score))
    
    # Almacenar p-value en la matriz
    p_values_matrix[i, j] <- p_value
  }
}


# # Calcular z-scores para toda la matriz
# all_z_scores <- as.vector((difference_values - mean(difference_values)) / sd(difference_values))
# 
# # Crear un marco de datos con los datos
# df_z_scores <- data.frame(z_scores = all_z_scores)
# 
# # Graficar la densidad de los z-scores
# ggplot(df_z_scores, aes(x = z_scores)) +
#   geom_density(fill = "blue", alpha = 0.5) +
#   labs(title = "Densidad de Z-scores", x = "Z-score", y = "Densidad")

# Convertir la matriz de p-values en un vector
all_p_values <- as.vector(p_values_matrix)

# Crear un marco de datos con los datos
df_p_values <- data.frame(p_values = all_p_values)

# Graficar la densidad de los p-values
ggplot(df_p_values, aes(x = p_values)) +
  geom_density(fill = "green", alpha = 0.5) +
  labs(title = "Densidad de P-values", x = "P-value", y = "Densidad")


# Imprimir los primeros resultados
#head(p_values_matrix)



# ggplot() +
#   geom_density(aes(x = as.vector(difference_values)), fill = "skyblue", color = "black") +
#   labs(title = " ", x = "cv", y = "Density")




# Crear un vector con los valores de cv_values (puedes ajustar esto segC:n tus necesidades)
#cv_vector <- as.vector(difference_values)


# 
# iniciales <- list(df = 3)
# 
# # Ajustar los datos simulados a una distribuciC3n chi-cuadrado
# ajuste_chi2 <- fitdist(cv_vector, "chisq", start = iniciales)
# 
# # Imprimir los resultados del ajuste
# print(ajuste_chi2)
# 
# # Crear un histograma de los datos simulados y la densidad ajustada
# hist(cv_vector, freq = FALSE, col = "lightblue", main = "Histograma y Ajuste a Chi-cuadrado", xlab = "Valor")
# curve(dchisq(x, df = 3), add = TRUE, col = "red", lwd = 2)
# 
# # Agregar leyenda
# legend("topright", legend = c("Datos Simulados", "Ajuste Chi-cuadrado"), col = c("lightblue", "red"), lwd = 2)
# 


#print(difference_values)
#load("./Data/results_data_Flevoland_100_5.Rdata")
#shapiro.test(as.vector(cv_values))
# Ejemplo para ajustar una distribuciC3n gamma
#fit_gamma <- fitdist(as.vector(cv_values), "gamma")

#summary(fit_gamma)

 # ggplot() +
 #   geom_density(aes(x = as.vector(difference_values)), fill = "skyblue", color = "black") +
 #   labs(title = " ", x = "cv", y = "Density")


 # # Calcular estadC-sticas descriptivas para la matriz cv_values
 # mean_val <- mean(as.vector(cv_values))
 # sd_val <- sd(as.vector(cv_values))
 # var_val <- var(as.vector(cv_values))
 # skewness_val <- e1071::skewness(as.vector(cv_values))
 # kurtosis_val <- e1071::kurtosis(as.vector(cv_values))
 # ad_test_result <- ad.test(as.vector(cv_values))
 # ad_p_value <- ad_test_result$p.value
 # cv_val <- sd_val / mean_val
 # 
 # # Imprimir o utilizar los valores calculados
 # cat("Media:", mean_val, "\n")
 # cat("DesviaciC3n EstC!ndar:", sd_val, "\n")
 # cat("Varianza:", var_val, "\n")
 # cat("AsimetrC-a:", skewness_val, "\n")
 # cat("Curtosis:", kurtosis_val, "\n")
 # cat("p-valor del test Anderson-Darling:", ad_p_value, "\n")
 # cat("Coeficiente de VariaciC3n:", cv_val, "\n")
 





# df <- data.frame(TestDifference = test_difference_vector)
# 
# 
# 
# ggplot(df, aes(x = TestDifference)) +
#   geom_density(color = "blue", fill = "lightblue", alpha = 0.7) +
#   labs(x = "Test Difference", y = "Density") +
#   theme_minimal() +
#   theme(
#     text = element_text(family = "serif"),
#     legend.position = "bottom"
#   )
# 
# 
# 
# # 
# mean_val <- mean(test_difference_vector)
# 
# sd_val <- sd(test_difference_vector)
# var_val <- var(test_difference_vector)
# skewness_val <- skewness(test_difference_vector)
# kurtosis_val <- kurtosis(test_difference_vector)
# ad_p_value <- ad.test(test_difference_vector)$p.value
# cv_val <- sd_val / mean_val
# 
# # #
# # cat("Mean:", mean_val, "\n")
# # 
# # cat("Standard Deviation:", sd_val, "\n")
# # cat("Variance:", var_val, "\n")
# # cat("Coefficient of Variation:", cv_val, "\n")
# # cat("Skewness:", skewness_val, "\n")
# # cat("Kurtosis:", kurtosis_val, "\n")
# # cat("Anderson-Darling p-value:", ad_p_value, "\n")
# 
# # 
# summary_stats <- data.frame(
#   Mean = mean_val,
#   SD = sd_val,
#   Variance = var_val,
#   CV = cv_val,
#   Skewness = skewness_val,
#   Kurtosis = kurtosis_val,
#   adpvalue = ad_p_value
# )
# 
# # qq_plot <- qqnorm(test_difference_vector, plot.it = FALSE)
# # qq_data <- data.frame(Theoretical = qq_plot$x, Sample = qq_plot$y)
# # 
# # 
# # ggplot(qq_data, aes(x = Theoretical, y = Sample)) +
# #   geom_point() +
# #   geom_abline(intercept = mean(test_difference_vector), slope = sd(test_difference_vector), col = "red", lty = 2) +
# #   labs(title = "QQ Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
# #   theme_minimal()
# 
# 
# save(summary_stats, file = "./Data/results.Rdata")
# 
# 
# calculate_p_value <- function(test_difference, sigma) {
#   epsilon <- (test_difference ) / sigma
#   p_value <- 2 * (1 - pnorm(abs(epsilon)))
#   
#   return(p_value)
# }
# 
# 
# 
# 
# #Calcular p-values para cada valor en test_difference_vector
# p_values <- sapply(test_difference_vector, function(value) {
#   calculate_p_value(value, sd_val )
# })
# 
# 
# p_values_df <- data.frame(
#   TestDifference = test_difference_vector,
#   PValue = p_values
# )
# 
# save(p_values,  file = "./Data/pvalues_Ottawa_100_5.Rdata")
# #print(p_values_df)
# 
# # ggplot(p_values_df, aes(x = PValue)) +
# #   geom_density(color = "blue", fill = "lightblue", alpha = 0.7) +
# #   labs(x = "Test Difference", y = "Density") +
# #   theme_minimal() +
# #   theme(
# #     text = element_text(family = "serif"),
# #     legend.position = "bottom"
# #   )
# 

