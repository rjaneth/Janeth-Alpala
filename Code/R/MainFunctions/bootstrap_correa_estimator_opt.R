bootstrap_correa_estimator_opt <- function(x, B) {
  v.Bootstrap <- numeric(B)
  
  for (b in seq_len(B)) {
    same_values <- TRUE
    while (same_values) {
      bx <- sample(x, replace = TRUE)
      if (!all(bx == bx[1])) {
        same_values <- FALSE
        entropy_result <- adaptive_correa_estimator(bx)
        if (is.finite(entropy_result)) {
          v.Bootstrap[b] <- entropy_result
          #cat("Entropy for Bootstrap Sample", b, ": ", v.Bootstrap[b], "\n")
        }
      }
    }
  }
  
  t <- adaptive_correa_estimator(x)
  estimated_mean <- mean(v.Bootstrap[!is.na(v.Bootstrap)])
  
  return(2 * t - estimated_mean)
}

# # Ejemplo de uso
# set.seed(123)  # Para reproducibilidad
# X <- rnorm(1000)  # Reemplaza esto con tu conjunto de datos
# B <- 1000
# system.time(result <- bootstrap_yee_estimator_optimized(X, B))
# print(paste("Resultado:", result))
