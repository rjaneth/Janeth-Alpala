adaptive_correa_estimator <- function(data) {
  n <- length(data)
  X_sorted <- sort(data)
  C_mn <- numeric(n-1)
  
  for (m in 1:(n-1)) {
    C_mn[m] <- 0
    
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
      
      C_mn[m] <- C_mn[m] - log(pmax(b_i, 0))
    }
    
    C_mn[m] <- C_mn[m] / n
  }
  
  estimated_entropy <- mean(C_mn)
  
  return(estimated_entropy)
}
