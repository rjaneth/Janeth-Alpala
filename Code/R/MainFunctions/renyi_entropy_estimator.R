renyi_entropy_estimator <- function(data, a) {
  n <- length(data)
  m <- round(sqrt(n) + 0.5)  # Espaciado m
  data_sorted <- sort(data)   # Ordenar los datos
  sum_term <- 0
  
  for (i in 1:n) {
    # Calcular el valor de ci
    if (i <= m) {
      ci <- 1 + (i - 1) / m
    } else if (i >= n - m + 1) {
      ci <- 1 + (n - i) / m
    } else {
      ci <- 2
    }
    
    # Manejar los límites para evitar salirse del rango
    if (i <= m) {
      diff_term <- data_sorted[i + m] - data_sorted[1]
    } else if (i >= n - m + 1) {
      diff_term <- data_sorted[n] - data_sorted[i - m]
    } else {
      diff_term <- data_sorted[i + m] - data_sorted[i - m]
    }
    
    # Sumar los términos al estimador usando la fórmula corregida
    sum_term <- sum_term + ((n / (ci * m)) * diff_term)^(1 - a)
  }
  
  # Calcular la entropía de Rényi usando la fórmula general
  renyi_entropy <- (1 / (1 - a)) * log(sum_term / n)
  
  return(renyi_entropy)
}




# renyi_entropy_estimator <- function(data, a) {
#   n <- length(data)
#   m <- round(sqrt(n) + 0.5)  # Espaciado m
#   data_sorted <- sort(data)   # Ordenar los datos
#   sum_term <- 0
#   
#   for (i in 1:n) {
#     # Calcular el valor de ci
#     if (i <= m) {
#       ci <- 1 + (i - 1) / m
#     } else if (i >= n - m + 1) {
#       ci <- 1 + (n - i) / m
#     } else {
#       ci <- 2
#     }
#     
#     # Calcular el término diferencial
#     if (i <= m) {
#       diff_term <- data_sorted[i + m] - data_sorted[1]
#     } else if (i >= n - m + 1) {
#       diff_term <- data_sorted[n] - data_sorted[i - m]
#     } else {
#       diff_term <- data_sorted[i + m] - data_sorted[i - m]
#     }
#     
#     # Sumar los términos al estimador
#     sum_term <- sum_term + ((n / (ci * m)) * (diff_term^(1 - a)))
#   }
#   
#   # Calcular la entropía de Rényi
#   renyi_entropy <- (1 / (1 - a)) * log(sum_term / n)
#   
#   return(renyi_entropy)
# }

#a_values <- c(0.01,0.1, 0.3, 0.5, 0.7, 0.8, 0.9,0.999,1.001, 1.1, 1.3, 1.5, 2, 3, 4, 5, 8)

# Ejemplo de uso
# data <- rnorm(100)  # Datos de ejemplo (puedes usar tus propios datos)
# renyi_entropies <- renyi_entropy_estimator(data, 0.999)
# renyi_entropies
