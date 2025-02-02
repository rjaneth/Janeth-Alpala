# enfoque numerico funciona
renyi_entropy_GI0_numeric <- function(alpha, mu, L, lambda) {
  if (lambda <= 0 || lambda == 1) {
    stop("Lambda must be greater than 0 and not equal to 1.")
  }
  
  # Definir gamma en funciC3n de mu y alpha
  gamma <- -mu * (alpha + 1)
  
  # Constante de la funciC3n de densidad
  C <- (L^L * gamma(-alpha + L)) / (gamma^alpha * gamma(-alpha) * gamma(L))
  
  # Definir la funciC3n f(z) elevada a la potencia lambda
  fz_lambda <- function(z) {
    fz <- C * z^(L - 1) / ( (gamma + L * z)^(L - alpha) )
    return(fz^lambda)
  }
  
  # Calcular la integral numC)ricamente
  integral_value <- integrate(fz_lambda, lower = 0, upper = Inf, subdivisions = 1000, rel.tol = .Machine$double.eps^0.25)
  
  if (integral_value$message != "OK") {
    warning("La integraciC3n numC)rica no fue exitosa: ", integral_value$message)
  }
  
  # Calcular la entropC-a de RC)nyi
  entropy <- (1 / (1 - lambda)) * log(integral_value$value)
  
  return(entropy)
}


# ParC!metros
alpha <- -3
mu <- 1.0
L <- 5
lambda <- 0.8

# Calcular la entropC-a numC)ricamente
entropy_value <- renyi_entropy_GI0_numeric(alpha, mu, L, lambda)
print(entropy_value)
