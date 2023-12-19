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
library(ggthemes)
library(ggthemes)
theme_set(theme_bw()  +
            theme(text=element_text(family="serif"),
                  legend.position = "top")
)


source("../MainFunctions/gamma_sar_sample.r")
source("../MainFunctions/entropy_gamma_sar.r")
source("../MainFunctions/entropy_gI0.r")
source("../MainFunctions/gi0_sample.r")

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
#The next function contains the functions: generate_samples, calculate_bias_mse, generate_plot
source("../Programs/functions_sample_bias_mse.R")# 

set.seed(1234567890, kind = "Mersenne-Twister")

sample_sizes <- c( 121)
R <- 1000
mu <- 100
L <- 5
alpha <- -10

# Function to generate samples for a given sample size and replication
generate_samples <- function(sample_size, replication, mu, alpha, L) {
  samples <- vector("list", replication)
  for (r in 1:replication) {
    samples[[r]] <- gi0_sample(mu, alpha, L, sample_size)
  }
  return(samples)
}

# Function to calculate Variance for both non-parametric and bootstrap estimators
calculate_entropy <- function(sample_sizes, R, mu, alpha, L) {
  estimators <- list(
    "Correa Bootstrap" = bootstrap_correa_estimator
  )
  
  output <- data.frame(SampleSize = integer(0), Estimator = character(0), Entropy = numeric(0))
  
  for (ssize in sample_sizes) {
    samples <- generate_samples(ssize, R, mu, alpha, L)
    
    for (estimator_name in names(estimators)) {
      estimator <- estimators[[estimator_name]]
      
      for (r in 1:R) {
        sample <- samples[[r]]
        
        if (grepl("Bootstrap", estimator_name)) {
          entropy_value <- estimator(sample, B = 10)
        } else {
          entropy_value <- estimator(sample)
        }
        
        output <- rbind(output, data.frame(SampleSize = ssize, Estimator = estimator_name, Entropy = entropy_value))
      }
    }
  }
  
  return(output)
}


empirical_distribution <- calculate_entropy(sample_sizes, R, mu, alpha, L)




ggplot(empirical_distribution, aes(x = Entropy, fill = Estimator)) +
  geom_histogram(position = "identity", alpha = 0.7, bins = 30) +
  facet_grid(. ~ SampleSize) +
  labs(title = "",
       x = "Entropy",
       y = "Frequency/Density") +
  theme_minimal()



