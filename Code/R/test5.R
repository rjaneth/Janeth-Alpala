rm(list = ls())

library(gtools)
library(stats4)
library(rmutil)
library(gsl)

mu <- 3
L <- 1
alphap <- -5

r <- 1000
#n_values <- c(9, 25, 49, 81, 200, 500, 1000)
n <- 1000

entropy_gamma_sar <- function(L, mu) {
  term1 <- L
  term2 <- log(gamma(L))
  term3 <- (1 - L) * digamma(L)
  term4 <- log(L / mu)
  
  entropy <- term1 +term2 + term3 - term4 
  
  return(entropy)
}


entropy_gI0 <- function(mu, alpha, L) {
  # entropía de GI0
  term1 <- -log(alphap / (mu * (1 + alphap))) - (1 - alphap) * digamma(-alphap)
  term2 <- log(-alphap / L) + (L - alphap) * digamma(L - alphap)
  term3 <- log(beta(L, -alphap))
  term4 <- (1 - L) * digamma(L)
  
  entropy <- term1 + term2 + term3 + term4
  return(entropy)
}

#muestras Gamma SAR
gamma_sar_sample <- function(L, mu, n) {
  samples <- rgamma(n, shape = L, rate = L / mu)
  return(samples)
}

#  estimador Van Es
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

D0 <- entropy_gamma_sar(L, mu)  

# ---
alpha <- 0.05
z_alpha_2 <- qnorm(1 - alpha / 2)
p_values <- c()
rechazo_h0 <- c()
means_entropies <- c()

entropies <- c()

  
  for (i in 1:r) {
    
    sample_data <- gamma_sar_sample(L, mu, n)
    entropy <- van_es_estimator(sample_data)
    entropies <- c(entropies, entropy)
    #log_var <- var(entropies)
    
    
    log_sample_data <- log(dgamma(sample_data, shape = L, rate = L / mu))
    log_var <- var(log_sample_data)
    
    
  }
  
  mean_entropy <- mean(entropies)
  means_entropies <- c(means_entropies, mean_entropy)
  #log_variances <- c(log_variances, log_var)

  
   #test statistics
    Z <- (sqrt(n) * (mean_entropy - D0)) / sqrt(log_var)
     
    p_value <- 2 *(1 - pnorm(abs(Z)))
    rechazar_h0 <- p_value < alpha
    
    hist(entropies, breaks = 30, probability = TRUE, col = "gray", border = "black",
         main = "",
         xlab = "Entropy", ylab = "Density")
    
    #hist(entropies, breaks = 30, probability = TRUE, col = "blue", main = "Non-Parametric Entropies Gamma Sar",
         #xlab = "Entropy", ylab = "Probability Density")
    grid()
                    
    # results
    cat("Sample Size:", n, "\n")
    cat("Z Statistic:", Z, "\n")
    cat("p_value:", p_value, "\n")
    cat("reject H0:", rechazar_h0, "\n\n")
    
    
    
    cat("media entropias:", mean_entropy, "\n")
    
    entropia_resultante <- entropy_gamma_sar(L, mu)
    cat("Entropía analítica Gamma:", entropia_resultante, "\n")
    
    
    entropia_resultante2 <- entropy_gI0(mu, alpha, L)
    cat("Entropía analítica GI0:", entropia_resultante2, "\n")
