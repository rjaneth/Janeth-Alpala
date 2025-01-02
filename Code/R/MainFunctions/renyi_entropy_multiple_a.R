# Función para calcular la entropía de Rényi para un conjunto de valores de 'a'
renyi_entropy_multiple_a <- function(data, a_values) {
  n <- length(data)
  m <- round(sqrt(n) + 0.5)  # Espaciado m
  data_sorted <- sort(data)   # Ordenar los datos
  entropies <- numeric(length(a_values))  # Vector para almacenar las entropías
  
  for (j in 1:length(a_values)) {
    a <- a_values[j]
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
      
      # Calcular el término diferencial
      if (i <= m) {
        diff_term <- data_sorted[i + m] - data_sorted[1]
      } else if (i >= n - m + 1) {
        diff_term <- data_sorted[n] - data_sorted[i - m]
      } else {
        diff_term <- data_sorted[i + m] - data_sorted[i - m]
      }
      
      # Sumar los términos al estimador
      sum_term <- sum_term + ((n / (ci * m)) * (diff_term^(1 - a)))
    }
    
    # Calcular la entropía de Rényi para el valor de a actual
    entropies[j] <- (1 / (1 - a)) * log(sum_term / n)
  }
  
  return(entropies)  # Devuelve un vector con las entropías para cada valor de a
}

# Conjunto de valores de 'a'
a_values <- c(0.01,0.1, 0.3, 0.5, 0.7, 0.8, 0.9,0.999,1.001, 1.1, 1.3, 1.5, 2, 3, 4, 5, 8)

# Ejemplo de uso
data <- rnorm(100)  # Datos de ejemplo (puedes usar tus propios datos)
renyi_entropies <- renyi_entropy_multiple_a(data, a_values)

# Mostrar los resultados
print(data.frame(a_values, renyi_entropies))
