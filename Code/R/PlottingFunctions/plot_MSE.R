rm(list = ls())
library(ggplot2)
library(invgamma)
library(ggthemes)
library(latex2exp)
source("../MainFunctions/gamma_sar_sample.r")
source("../MainFunctions/entropy_gamma_sar.r")

source("../MainFunctions/van_es_estimator.r")
source("../MainFunctions/correa_estimator.r")
source("../MainFunctions/ebrahimi_estimator.r")
source("../MainFunctions/noughabi_arghami_estimator.r")
source("../MainFunctions/vasicek_estimator.r")
source("../MainFunctions/al_omari_1_estimator.r")
source("../MainFunctions/al_omari_2_estimator.r")


sample_sizes <- c(9, 25, 49, 81, 121, 225)
R <- 1000
mu <- 1
L <- 2

# True entropy value for the gamma SAR distribution
true_entropy <- entropy_gamma_sar(L, mu)



output_mse <- data.frame(SampleSize = integer(0), Estimator = character(0), MSE = numeric(0))
#List of estimator
estimators <- list(
  list(func = van_es_estimator, label = "Van Es"),
  list(func = correa_estimator, label = "Correa"),
  list(func = ebrahimi_estimator, label = "Ebrahimi"),
  list(func = noughabi_arghami_estimator, label = "Noughabi Arghami"),
  list(func = vasicek_estimator, label = "Vasicek"),
  list(func = al_omari_1_estimator, label = "Al Omari 1"),
  list(func = al_omari_2_estimator, label = "Al Omari 2")
)



for (estimator in estimators) {
  for (ssize in sample_sizes) {
    squared_errors <- numeric(R)
    for (r in 1:R) {
      sample <- gamma_sar_sample(L, mu, ssize)
      estimated_entropy <- estimator$func(sample)
      squared_errors[r] <- (estimated_entropy - true_entropy)^2
    }
    
    # Calculate the MSE
    mse <- mean(squared_errors)
    
    output_mse <- rbind(output_mse, data.frame(SampleSize = ssize, Estimator = estimator$label, MSE = mse))
  }
}

# Plot the MSE as a function of the sample size
ggplot(output_mse, aes(x = SampleSize, y = MSE, color = Estimator)) +
  geom_line() +
  geom_point() +
  xlab("Sample Size") +
  ylab("MSE") 

  
