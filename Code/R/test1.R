
rm(list = ls())
library(gtools)
library(stats4)
library(rmutil)
library(gsl)

# mu <-100
# L <- 18
# alpha <- -100


mu <- 5
L <- 1
alphap <- -4

r <- 10000
n_values <- c(1000)  


# Function 
entropy_gamma_sar<- function(L, mu) {
  term1 <- L
  term2 <- log(gamma(L))
  term3 <- (1 - L) * digamma(L)
  term4 <- log(L / mu)
  
  entropy <- term1 +term2 + term3 - term4 
  
  return(entropy)
}

entropy_gI0 <- function(mu, alphap, L) {
  # entropía de GI0
  term1 <- -log(alphap / (mu * (1 + alphap))) - (1 - alphap) * digamma(-alphap)
  term2 <- log(-alphap / L) + (L - alphap) * digamma(L - alphap)
  term3 <- log(beta(L, -alphap))
  term4 <- (1 - L) * digamma(L)
  
  entropy <- term1 + term2 + term3 + term4
  return(entropy)
}


gamma_sar_sample <- function(L, mu, n) {
  samples <- rgamma(n, shape = L, rate = L / mu)
  return(samples)
}



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


means_entropies <- c()

for (n in n_values) {
  entropies <- c()
  for (i in 1:r) {
    sample_data <- gamma_sar_sample(L, mu, n)
    entropy <- van_es_estimator(sample_data)
    entropies <- c(entropies, entropy)
    log_gamma_sar<- log(dgamma(sample_data, shape = L, rate = L / mu))
    #log_var <- var(log_sample_data)
    
  }
  
  mean_entropy <- mean(entropies)
  means_entropies <- c(means_entropies, mean_entropy)
}


# plot(n_values, means_entropies, type = "o", main = "Convergence of Mean Non-Parametric Entropies",
#      xlab = "Sample Size (n)", ylab = "Mean Entropy")
# grid()


entropy_Gamma_sar <- entropy_gamma_sar(L, mu)
cat("mean of Non-Parametric Entropies gamma sar :", means_entropies, "\n")
cat("Analitical Entropy gamma sar:", entropy_Gamma_sar, "\n")

entropy_GI0 <- entropy_gI0(mu, alphap, L)
cat("Analitical Entropy GIO:", entropy_GI0, "\n")



alpha <- 0.05  
D0 <- entropy_Gamma_sar  


Z <- (means_entropies - D0) / (sqrt(var(log_gamma_sar)) / sqrt(n_values))

# Valor crítico
z_alpha_2 <- qnorm(1 - alpha/2)




p_values <- 2 * (1 - pnorm(abs(Z)))


rechazo <- p_values < alpha

# Resultados
resultados <- data.frame(
  Sample_Size = n_values,
  Mean_Entropy = means_entropies,
  Z_Statistic = Z,
  P_Value = p_values,
  Rechazo_H0 = rechazo
)

print(resultados)

