rm(list = ls())
# ctrl + L  Clear the Console
#setwd("C:/Users/luiso/Documents/Github/Janeth-Alpala/Code/R")
library(gtools)
library(stats4)
library(rmutil)
library(gsl)

# mu <-100
# L <- 18
# alpha <- -100


mu <- 3
L <- 2
alpha <- -5

r <- 10000
n_values <- c(9, 25, 49, 81, 121, 500, 1000, 5000, 10000 )  


# Function 
entropy_gamma_sar<- function(L, mu) {
  term1 <- L
  term2 <- log(gamma(L))
  term3 <- (1 - L) * digamma(L)
  term4 <- log(L / mu)
  
  entropy <- term1 +term2 + term3 - term4 
  
  return(entropy)
}

entropy_gI0 <- function(mu, alpha, L) {
  # entropía de GI0
  term1 <- -log(alpha / (mu * (1 + alpha))) - (1 - alpha) * digamma(-alpha)
  term2 <- log(-alpha / L) + (L - alpha) * digamma(L - alpha)
  term3 <- log(beta(L, -alpha))
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
    
  }
  
  mean_entropy <- mean(entropies)
  means_entropies <- c(means_entropies, mean_entropy)
}

#Convergence of Mean Non-Parametric Entropies
plot(n_values, means_entropies, type = "o", main = "",
     xlab = "Sample Size (n)", ylab = "Mean Entropy")
grid()

# hist(entropies, breaks = 30, probability = TRUE, col = "blue", main = "Histogram of Non-Parametric Entropies Gamma Sar",
#      xlab = "Entropy", ylab = "Probability Density")
# grid()
# 
entropy_Gamma_sar <- entropy_gamma_sar(L, mu)
cat("mean of Non-Parametric Entropies gamma sar :", means_entropies, "\n")
cat("Analitical Entropy gamma sar:", entropy_Gamma_sar, "\n")

entropy_GI0 <- entropy_gI0(mu, alpha, L)
cat("Analitical Entropy GIO:", entropy_GI0, "\n")

#cat("vector:", entropies, "\n")
#cat("vector2:", entropy, "\n")
#cat("vector3:", sample_data, "\n")
