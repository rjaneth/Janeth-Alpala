# Función para calcular la entropía de Rényi de la distribución G_I^0 con sustitución de gamma
entropy_renyi_GI0_2 <- function(lambda, alpha, mu, L) {
  # Validación básica de los parámetros
  if (lambda <= 0 || lambda == 1) {
    stop("El parámetro 'lambda' debe ser mayor que 0 y diferente de 1.")
  }
  if (L <= 0) {
    stop("El parámetro 'L' debe ser mayor que 0.")
  }
  if (mu <= 0) {
    stop("El parámetro 'mu' debe ser mayor que 0.")
  }
  if (alpha >= -1) {
    stop("El parámetro 'alpha' debe ser menor que -1.")
  }
  
  # Cálculo de gamma utilizando la sustitución gamma = -mu * (alpha + 1)
  gamma <- -mu * (alpha + 1)
  
  # Verificar que gamma sea positivo
  if (gamma <= 0) {
    stop("El valor calculado de 'gamma' debe ser mayor que 0. Verifique que 'mu' > 0 y 'alpha' < -1.")
  }
  
  # Cálculo de a y b
  a <- lambda * (L - 1) + 1
  b <- lambda * (1 - alpha) - 1
  
  # Cálculo de los términos
  term1 <- lambda * lgamma(L - alpha)
  term2 <- - lambda * lgamma(-alpha)
  term3 <- - lambda * lgamma(L)
  
  # Reemplazar ln(gamma) por ln(mu) + ln(- (alpha + 1))
  ln_gamma <- log(mu) + log(- (alpha + 1))
  term4 <- (lambda - 1) * (ln_gamma - log(L))
  
  term5 <- lgamma(a)
  term6 <- lgamma(b)
  term7 <- - lgamma(lambda * (L - alpha))
  
  # Cálculo de la entropía de Rényi
  H <- (term1 + term2 + term3 + term4 + term5 + term6 + term7) / (1 - lambda)
  
  return(H)
}


# Parámetros
lambda <- 0.8
alpha <- -3
mu <- 1
L <- 5

# Cálculo de la entropía de Rényi
H <- entropy_renyi_GI0_2(lambda, alpha, mu, L)
print(H)


