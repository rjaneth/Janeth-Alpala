al_omari_estimator_m <- function(X) {
  n <- length(X)
  H_AO <- numeric(n-1)
  
  sorted_X <- sort(X)
  
  for (m in 1:(n-1)) {
    sum_term <- 0
    
    for (i in 1:n) {
      if (i <= m) {
        omega_i <- 3/2
        diff_term <- sorted_X[i + m] - sorted_X[1]
      } else if (i >= n - m + 1) {
        omega_i <- 3/2
        diff_term <- sorted_X[n] - sorted_X[i - m]
      } else {
        omega_i <- 2
        diff_term <- sorted_X[i + m] - sorted_X[i - m]
      }
      
      # Agregar comprobación para evitar logaritmo de cero y manejar NA
      if (!is.na(diff_term) && diff_term != 0) {
        sum_term <- sum_term + log((n / (omega_i * m)) * diff_term)
      }
    }
    
    H_AO[m] <- sum_term / n
  }
  
  estimated_entropy <- (1/(n-1)) * sum(H_AO)
  
  return(estimated_entropy)
}

# set.seed(123)
# test_data <- rnorm(1000)
# result <- al_omari_estimator(test_data)
# print(result)
