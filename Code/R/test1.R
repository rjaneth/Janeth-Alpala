
rm(list = ls())
library(gtools)
library(stats4)
library(rmutil)
library(gsl)
library(invgamma)



alpha <- -2
mu <- 1
L <- 1

gam=-mu*(alpha+1)

r <- 10000
n_values <- c(9, 25, 49, 81, 200,1000, 5000,  10000) 


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


means_entropies <- c()
means_entropies2 <- c()
entropies2 <- c()
for (n in n_values) {
  entropies <- c()
  entropies2 <- c()
  for (i in 1:r) {
    sample_data <- gamma_sar_sample(L, mu, n)
    GI0_samples <- Z_samples(L, alpha, mu,  n)
    Np_entropy_gamma <- van_es_estimator(sample_data)
    NP_entropy_GI0 <- van_es_estimator(GI0_samples)
    entropies <- c(entropies, Np_entropy_gamma)
    entropies2 <- c(entropies2, NP_entropy_GI0)
    log_gamma_sar<- log(dgamma(sample_data, shape = L, rate = L / mu))
    #log_var <- var(log_sample_data)
    
  }
  #Gamma SAR
  mean_entropy <- mean(entropies)
  means_entropies <- c(means_entropies, mean_entropy)
  # GIO
  mean_entropy2 <- mean(entropies2)
  means_entropies2 <- c(means_entropies2, mean_entropy2)
}


plot(n_values, means_entropies, type = "o", main = "Convergence of Mean Non-Parametric Entropies",
     xlab = "Sample Size (n)", ylab = "Mean Entropy")
grid()


entropy_Gamma_sar <- entropy_gamma_sar(L, mu)
cat("mean of Non-Parametric Entropies gamma sar :", means_entropies, "\n")
cat("Analitical Entropy gamma sar:", entropy_Gamma_sar, "\n")

entropy_GI0 <- entropy_gI0(mu, alpha, L)
cat("Analitical Entropy GIO:", entropy_GI0, "\n")

cat("mean of Non-Parametric Entropies GI0 :", means_entropies2, "\n")

second <- second_term(alpha, L)
cat("second term of GI0 :", second, "\n")



# # test
#   alpha1 <- 0.05
#   D0 <- entropy_Gamma_sar
# 
# 
#   Z <- (means_entropies - D0) / (sqrt(var(log_gamma_sar)) / sqrt(n_values))
# 
#   # Valor crítico
#   z_alpha_2 <- qnorm(1 - alpha1/2)
# 
# 
# 
# 
#   p_values <- 2 * (1 - pnorm(abs(Z)))
# 
# 
#   rechazo <- p_values < alpha1
# 
#   # Resultados
#   resultados <- data.frame(
#     Sample_Size = n_values,
#     Mean_Entropy = means_entropies,
#     Z_Statistic = Z,
#     P_Value = p_values,
#     Rechazo_H0 = rechazo
#   )
# 
#   print(resultados)
# 
