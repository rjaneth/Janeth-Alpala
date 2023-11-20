rm(list = ls())
library(ggplot2)
library(tidyr)
library(gtools)
library(stats4)
library(rmutil)
library(invgamma)
library(ggplot2)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ggsci)
library(wesanderson)
library(ggpubr)

# Load the source files (make sure you have the correct file paths)
source("../MainFunctions/gamma_sar_sample.r")
source("../MainFunctions/entropy_gamma_sar.r")
source("../MainFunctions/van_es_estimator.r")
source("../MainFunctions/correa_estimator.r")
source("../MainFunctions/ebrahimi_estimator.r")
source("../MainFunctions/noughabi_arghami_estimator.r")
source("../MainFunctions/vasicek_estimator.r")
source("../MainFunctions/al_omari_1_estimator.r")
source("../MainFunctions/al_omari_2_estimator.r")
source("../MainFunctions/bootstrap_van_es_estimator.r")
source("../MainFunctions/bootstrap_correa_estimator.r")
source("../MainFunctions/bootstrap_ebrahimi_estimator.r")
source("../MainFunctions/bootstrap_noughabi_arghami_estimator.r")
source("../MainFunctions/bootstrap_vasicek_estimator.r")
source("../MainFunctions/bootstrap_al_omari_1_estimator.r")
source("../MainFunctions/bootstrap_al_omari_2_estimator.r")

set.seed(1234567890, kind = "Mersenne-Twister")

sample_sizes <- c( 9, 25, 49, 81, 121, 225)
R <- 100
mu <- 1
L <- 1

# Function to calculate Bias and MSE for both non-parametric and bootstrap estimators
calculate_bias_mse <- function(sample_sizes, R, mu, L) {
  true_entropy <- entropy_gamma_sar(L, mu)
  
  # Define a list of estimators
  estimators <- list(
    "Van Es" = van_es_estimator,
    "Correa" = correa_estimator,
    "Ebrahimi" = ebrahimi_estimator,
    "Noughabi Arghami" = noughabi_arghami_estimator,
    "Vasicek" = vasicek_estimator,
    "Al Omari 1" = al_omari_1_estimator,
    "Al Omari 2" = al_omari_2_estimator,
    "Van Es Bootstrap" = bootstrap_van_es_estimator,
    "Correa Bootstrap" = bootstrap_correa_estimator,
    "Ebrahimi Bootstrap" = bootstrap_ebrahimi_estimator,
    "Noughabi Arghami Bootstrap" = bootstrap_noughabi_arghami_estimator,
    "Vasicek Bootstrap" = bootstrap_vasicek_estimator,
    "Al Omari 1 Bootstrap" = bootstrap_al_omari_1_estimator,
    "Al Omari 2 Bootstrap" = bootstrap_al_omari_2_estimator
  )
  
  output <- NULL
  
  for(ssize in sample_sizes) {
    for(estimator_name in names(estimators)) {
      estimator <- estimators[[estimator_name]]
      
      v.entropy <- NULL
      
      for(r in 1:R) {
        sample <- gamma_sar_sample(L, mu, ssize)
        if (grepl("Bootstrap", estimator_name)) {
          v.entropy[r] <- estimator(sample, B = 10)
        } else {
          v.entropy[r] <- estimator(sample)
        }
      }
      
      mse <- mean((v.entropy - true_entropy)^2)
      bias <- mean(v.entropy) - true_entropy
      
      output <- rbind(output, data.frame(SampleSize = ssize, Estimator = estimator_name, Bias = bias, MSE = mse))
    }
  }
  
  return(output)
}



results <- calculate_bias_mse(sample_sizes, R, mu, L)

# Print the results
print(results)

df <- as.data.frame(results)

# Plot Bias and MSE in separate subplots
plot_bias <- ggplot(df, aes(x = SampleSize, y = Bias, color = Estimator)) +
  geom_hline(yintercept = 0) +
  geom_point(size = 2) +
  geom_line(linetype = "solid", size = 0.5) +
  labs(title = "",
       y = "Bias", x = "Sample size") +
  guides(color = guide_legend(title = "Estimator"))

plot_mse <- ggplot(df, aes(x = SampleSize, y = MSE, color = Estimator)) +
  geom_hline(yintercept = 0) +
  geom_point(size = 2) +
  geom_line(linetype = "solid", size = 0.5) +
  labs(title = "",
       y = "MSE", x = "Sample size") +
  guides(color = guide_legend(title = "Estimator"))

# Arrange both plots in a single figure with subplots
ggarrange(plot_bias, plot_mse, ncol=2, nrow=1, common.legend = TRUE, legend="top")
#grid.arrange(plot_bias, plot_mse, ncol = 2)
