
rm(list = ls())
# ctrl + L  Clear the Console
#setwd("C:/Users/luiso/Documents/Github/Janeth-Alpala/Code/R")
library(gtools)
library(stats4)
library(rmutil)
library(gsl)



entropy_g0 <- function(alpha, gam, L) {
  # Calcular la entropía de G^0
  term1 <- -log(-alpha / gam) - (1 - alpha) * digamma(-alpha)
  term2 <- log(-alpha / L) + (L - alpha) * digamma(L - alpha)
  term3 <- log(beta(L, -alpha))
  term4 <- (1 - L) * digamma(L)
  
  entropy <- term1 + term2 + term3 + term4
  return(entropy)
}

rGI0= function(alpha,gam,n,L = 1,gen){
  switch(gen,
         "Pareto" = rpareto(n,-alpha,gam),
         "Gamma" = rgamma(n,L,L)/rgamma(n,-alpha,gam),
         "Chi2" = gam*rchisq(n,2)/rchisq(n,-2*alpha),
         "F-Snedecor" = (-1)*gam*rf(n,2,-2*alpha)/ alpha)
}




#  non-parametric entropy - Van Es estimator
van_es_estimator <- function(data) {
  n <- length(data)
  m <- floor(sqrt(n) + 0.5)  # m-spacing
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


alpha <- -2
gam <- 2
n <- 20000
gen <- "Pareto"  


L <- 1

# Number of samples 
r <- 10000


entropies <- vector()
for (i in 1:r) {
  sample <- rGI0(alpha, gam, n, L , gen)
  entropy <- van_es_estimator(sample)
  entropies <- c(entropies, entropy)
}

mean_entropy <- mean(entropies)


hist(entropies, breaks = 30, probability = TRUE, col = "blue", main = "Histogram of Non-Parametric Entropies G0",
     xlab = "Entropy", ylab = "Probability Density")
grid()

entropia_g0 <- entropy_g0(alpha, gam, L)
cat("Entropía de G^0:", entropia_g0, "\n")
cat("means_entropies:", mean_entropy, "\n")
