#rm(list = ls())
#load("./Data/results_Chicago_1024_9_AO_200b_36L.Rdata")
#save(difference_values, x, file = "./Data/results_Frankfurt_1024_9W_AO_100b_5L.Rdata")
#load("./Data/results_Frankfurt_1024_9W_AO_100b_5L.Rdata")
#load("./Data/results_Chicago_1024_9_AO_100b_6L.Rdata")
#save(difference_values, x, file = "./Data/results_mexico_600_7W_AO_50b_6L.Rdata")save(difference_values, x, file = "./Data/results_mexico_600_9W_AO_100b_6L.Rdata")
#load("./Data/results_mexico_600_9W_AO_100b_6L.Rdata")
#load("./Data/results_chicago_urban_1024_9W_AO_100b_6L.Rdata")
#load("./Data/results_Chicago_400_9W_AO_100b_5L.Rdata")
#load("./Data/results_Lake_512_9W_AO_100b_36L.Rdata")
#load("./Data/results_Chicago_1024_9_AO_200b_36L.Rdata")
#save(difference_values, x, file = "./Data/results_Chicago_400_9W_AO_100b_5L.Rdata")
#save(difference_values, x, file = "./Data/results_Lake_512_9W_AO_100b_36L.Rdata")
#load("./Data/results_ottawa_city_512_7_AO_300b.Rdata")
#load("./Data/results_data_simulated_z_up_le_entropy_250_7.Rdata")
#load("../Programs/Data/results_Panama_700_AO_7W_L5_200b.Rdata")
load("./Data/results_SanFranciscoN_400_renyi_b200_07_L1_w7.Rdata")
#load("./Data/results_Phantom_4_renyi_B100_w7_09_L5.Rdata")

calculate_p_values_and_epsilon <- function(data_matrix) {
  rows <- nrow(data_matrix)
  cols <- ncol(data_matrix)
  
  mu <- mean(data_matrix, na.rm = TRUE)  # Calcular la media
  sigma <- sd(data_matrix, na.rm = TRUE)  # Calcular la desviación estándar
  
  p_values_matrix <- matrix(NA, nrow = rows, ncol = cols)
  epsilon_matrix <- matrix(NA, nrow = rows, ncol = cols)
  
  for (i in 1:rows) {
    for (j in 1:cols) {
      test_difference <- data_matrix[i, j]
      
      epsilon <- (test_difference - 0) / sigma # le coloqué mu =0 y fucniona!
      epsilon_matrix[i, j] <- epsilon  # Guardar el valor de epsilon
      
      p_value <- 2 * pnorm(-abs(epsilon))
      
      p_values_matrix[i, j] <- p_value
    }
  }
  
  return(list(p_values_matrix = p_values_matrix, epsilon_matrix = epsilon_matrix))
}

# Llamada a la función y guardado de los resultados
result <- calculate_p_values_and_epsilon(difference_values)
p_values_matrix <- result$p_values_matrix
epsilon_matrix <- result$epsilon_matrix

# Guardar los resultados
#save(p_values_matrix, epsilon_matrix, file = "./Data/results_pvalue_epsilon.Rdata")

# Definir el umbral de significancia
significance_level <- 0.05

# Crear una matriz binaria de p-valores significativos
significant_pixels <- p_values_matrix <= significance_level

# Visualizar la matriz binaria
image(significant_pixels, col = c( "white","black"))  # Blanco para no significativo, negro para significativo

# Definir la escala de grises para el degradado
# n_shades <- 100  # Número de tonos de gris en el degradado
# gray_palette <- gray.colors(n_shades)  # Paleta de colores en escala de grises
# 
# # Convertir p-values a tonos de gris en función de su significancia
# gray_scale_values <- round((p_values_matrix / significance_level) * n_shades)  # Escalar los p-values a la escala de grises
# 
# # Limitar los valores escalados para evitar valores fuera del rango de la paleta de colores
# gray_scale_values[gray_scale_values > n_shades] <- n_shades
# gray_scale_values[gray_scale_values < 1] <- 1
# 
# # Visualizar la matriz de tonos de gris
#image(gray_scale_values, col = gray_palette)

