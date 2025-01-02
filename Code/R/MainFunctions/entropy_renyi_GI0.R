# Función para calcular la entropía de Rényi H_lambda(Z) para la distribución G_I^0
entropy_renyi_GI0 <- function(lambda, L, mu, alpha) {
  # Verificación de los parámetros de entrada
  if (lambda <= 0 || lambda == 1) {
    stop("El parámetro 'lambda' debe ser positivo y diferente de 1.")
  }
  if (L < 1) {
    stop("El parámetro 'L' debe ser mayor o igual a 1.")
  }
  if (mu <= 0) {
    stop("El parámetro 'mu' debe ser positivo.")
  }
  if (alpha >= -1) {
    stop("El parámetro 'alpha' debe ser menor que -1.")
  }
  if (- (alpha + 1) <= 0) {
    stop("El valor de '-(alpha + 1)' debe ser positivo para el logaritmo.")
  }
  
  # Verificación de los argumentos de las funciones gamma
  if (L - alpha <= 0) {
    stop("El valor de 'L - alpha' debe ser positivo para la función gamma.")
  }
  if (-alpha <= 0) {
    stop("El valor de '-alpha' debe ser positivo para la función gamma.")
  }
  if (lambda * (L - alpha) <= 0) {
    stop("El valor de 'lambda * (L - alpha)' debe ser positivo para la función gamma.")
  }
  if (-lambda * (alpha - 1) - 1 <= 0) {
    stop("El argumento de la función gamma debe ser positivo.")
  }
  if (lambda * (L - 1) + 1 <= 0) {
    stop("El valor de 'lambda * (L - 1) + 1' debe ser positivo para la función gamma.")
  }
  
  # Cálculo de los términos de la entropía
  term1 <- ((1 - lambda * L) / (1 - lambda)) * (log(mu) + log(- (alpha + 1)))
  term2 <- -log(L)
  term3 <- (lambda / (1 - lambda)) * (lgamma(L - alpha) - lgamma(-alpha) - lgamma(L))
  term4 <- (1 / (1 - lambda)) * (lgamma(lambda * (L - 1) + 1) +
                                   lgamma(-lambda * (alpha - 1) - 1) -
                                   lgamma(lambda * (L - alpha)))
  
  # Suma de todos los términos para obtener la entropía
  H_lambda_Z <- term1 + term2 + term3 + term4
  
  # Devolver el valor calculado de la entropía
  return(H_lambda_Z)
}

lambda <- 0.8
alpha <- -3
mu <- 1
L <- 5

# Cálculo de la entropía de Rényi
H <- entropy_renyi_GI0(lambda, L, mu, alpha)
print(H)