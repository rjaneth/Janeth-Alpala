bootstrap_renyi_entropy_estimator <- function(x, B, a) {
  
  v.Bootstrap <- rep(0, B)
  
  for(b in 1:B) {
    same_values <- TRUE
    while (same_values) {
      # Muestra bootstrap con reemplazo
      bx <- sample(x, replace = TRUE)
      if (!all(bx == bx[1])) {
        same_values <- FALSE
        entropy_result <- renyi_entropy_estimator(bx, a)
        if (is.finite(entropy_result)) {
          v.Bootstrap[b] <- entropy_result
        }
      }
    }
  }
  
  # Calcular la entropía de Rényi para la muestra original
  t <- renyi_entropy_estimator(x, a)
  
  # Calcular el promedio de los valores bootstrap
  estimated_mean <- mean(v.Bootstrap[!is.na(v.Bootstrap)])
  
  # Estimador mejorado
  return(2 * t - estimated_mean)
}
