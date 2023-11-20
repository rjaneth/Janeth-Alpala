rm(list = ls())

#options(digits = 5)

# Resto de tu código aquí

# Restaurar la opción digits a su valor original si es necesario
#on.exit(options(digits = getOption("digits")))
library(ggplot2)
library(tidyr)
library(gtools)
library(stats4)
library(rmutil)
library(invgamma)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ggsci)
library(wesanderson)
library(ggpubr)
library(patchwork)
theme_set(theme_bw()  +
            theme(text=element_text(family="serif"),
                  legend.position = "top", panel.grid = element_blank())# Gridtop , right , bottom , or left#, panel.grid = element_blank()
)


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

#The next function contains the functions: generate_samples, calculate_bias_mse, generate_plot
source("../Programs/functions_sample_bias_mse.R")# 

set.seed(1234567890, kind = "Mersenne-Twister")


# Define sample sizes for analysis
sample_sizes <- c(9, 25, 49, 81, 121, 225)

# Number of replications
R <- 10

# Number of bootstrap replications
B <- 10

# mu values 
#mu_values <- c(1)

 mu_values <- c(1, 2, 3, 4)

#L values
#mu <- 1
L <- 1

# Define a list of estimators
estimators <- list(
  "Van Es" = van_es_estimator,
  "Van Es Bootstrap" = bootstrap_van_es_estimator,
  "Correa" = correa_estimator,
  #"Ebrahimi" = ebrahimi_estimator,
  #"Noughabi Arghami" = noughabi_arghami_estimator,
  "Vasicek" = vasicek_estimator,
  #"Al Omari 1" = al_omari_1_estimator,
  #"Al Omari 2" = al_omari_2_estimator,
  "Correa Bootstrap" = bootstrap_correa_estimator,
  #"Ebrahimi Bootstrap" = bootstrap_ebrahimi_estimator,
  #"Noughabi Arghami Bootstrap" = bootstrap_noughabi_arghami_estimator,
  "Vasicek Bootstrap" = bootstrap_vasicek_estimator
  # "Al Omari 1 Bootstrap" = bootstrap_al_omari_1_estimator,
  #"Al Omari 2 Bootstrap" = bootstrap_al_omari_2_estimator
)

# Generate and plot the results
combined_plot <- generate_plot(sample_sizes, R, B, mu_values, L, estimators, ncol = 2, nrow = 2)

# Print the combined plot
print(combined_plot)