# Functions
# Entropy Gamma SAR
entropy_gamma_sar <- function(L, mu) {
  
  return(L - log(L) + log(gamma(L)) + (1 - L) * digamma(L) + log(mu))
  
}

# Entropy GI0
entropy_gI0 <- function(mu, alpha, L) {
  
  term1 <- L - log(L) + log(gamma(L)) + (1 - L) * digamma(L) + log(mu)   
  term2 <- -L - log(gamma(L-alpha)) + (L-alpha)*(digamma(L - alpha))- (1-alpha)*digamma(- alpha)+log(-1 - alpha)+log(gamma(-alpha))
  
  entropy <- term1 + term2 
  return(entropy)
}


# only second term of GI0
second_term <- function(alpha, L) {
  return(-L - log(gamma(L-alpha)) + (L-alpha)*(digamma(L - alpha))- (1-alpha)*digamma(- alpha)+log(-1 - alpha)+log(gamma(-alpha)))
  
}


#Samples Gamma SAR
gamma_sar_sample <- function(L, mu, n) {
  samples <- rgamma(n, shape = L, rate = L / mu)
  return(samples)
}


#Samples GI0
Z_samples <- function(L, alpha, mu,  n) {
  
  X_samples <- rinvgamma(n, shape = -alpha, rate =mu*(-alpha-1))
  
  Y_samples <- rgamma(n, shape = L, rate = L)
  Z_samples <- X_samples * Y_samples
  return(Z_samples)
}



# nonparametric estimator



van_es_estimator <- function(data) {
  n <- length(data)
  m <- round(sqrt(n) + 0.5)  # m-spacing
  data_sorted <- sort(data)
  sum_term1 <- 0
  sum_term2 <- 0
  
  for (i in 1:(n - m)) {
    sum_term1 <- sum_term1 + log(((n + 1) / m) * (data_sorted[i + m] - data_sorted[i]))
  }
  
  for (k in m:n) {
    sum_term2 <- sum_term2 + 1 / k
  }
  
  sum_term3 <- log(m / (n + 1))
  
  return((sum_term1 / (n - m)) + sum_term2+ sum_term3)
}
