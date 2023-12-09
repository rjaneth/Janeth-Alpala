correa_estimator_opt <- function(data) {
  n <- length(data)
  m <- round(sqrt(n) + 0.5)
  X_sorted <- sort(data)
  C_mn <- 0
  
  for (i in 1:n) {
    i_minus_m <- max(1, i - m)
    i_plus_m <- min(n, i + m)
    
    X_i <- X_sorted[i]
    X_window <- X_sorted[i_minus_m:i_plus_m]
    X_bar_i <- mean(X_window)
    
    differences <- X_window - X_bar_i
    num <- sum(differences * (i_minus_m:i_plus_m - i))
    den <- n * sum(differences^2)
    
    b_i <- ifelse(den == 0, 0, num / den)
    
    C_mn <- C_mn - log(pmax(b_i, 0))
  }
  
  C_mn <- C_mn / n
  
  return(C_mn)
}

# set.seed(123)  
# X <- rnorm(100000)  
# #system.time(entropy_estimate_large <- yee_entropy_estimator_optimized(X_large))
# system.time(entropy_estimate <- correa_estimator_3(X))
# print(paste("Estimación de entropía:", entropy_estimate))