

yee_estimator<- function(X) {
  n <- length(X)
  H_Y <- numeric(n-1)
  
  sorted_X <- sort(X)
  
  for (m in 1:(n-1)) {
    differences <- sorted_X[(m+1):n] - sorted_X[1:(n-m)]
    term_sum <- sum(log((n/m) * differences))
    
    H_Y[m] <- (1/(n-m)) * term_sum - digamma(m) + log(m)
  }
  
  estimated_entropy <- (1/(n-1)) * sum(H_Y)
  
  return(estimated_entropy)
}

























# yee_estimator<- function(X) {
#   n <- length(X)
#   H_Y <- numeric(n-1)
#   
#   sorted_X <- sort(X)
#   
#   for (m in 1:(n-1)) {
#     differences <- sorted_X[(m+1):n] - sorted_X[1:(n-m)]
#     term_sum <- sum(log((n/m) * differences))
#     
#     H_Y[m] <- (1/(n-m)) * term_sum - digamma(m) + log(m)
#   }
#   
#   estimated_entropy <- (1/(n-1)) * sum(H_Y)
#   
#   return(estimated_entropy)
# }


# set.seed(123)  # Para reproducibilidad
# X_large <- rnorm(1000)  # Reemplaza esto con tu conjunto de datos
# system.time(entropy_estimate_large <- yee_entropy_estimator_optimized(X_large))
# print(paste("Estimación de entropía:", entropy_estimate_large))
