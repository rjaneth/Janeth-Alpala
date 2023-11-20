
rm(list = ls())
library(ggplot2)
library(gtools)
library(stats4)
library(rmutil)
#library(gsl)
library(invgamma)
library(ggplot2)
library(ggthemes)
theme_set(theme_pander() +
            theme(text=element_text(family="serif"),
                  legend.position = "top")
)
source("../MainFunctions/gamma_sar_sample.r")
source("../MainFunctions/entropy_gamma_sar.r")

source("../MainFunctions/van_es_estimator.r")
source("../MainFunctions/correa_estimator.r")
source("../MainFunctions/ebrahimi_estimator.r")
source("../MainFunctions/noughabi_arghami_estimator.r")
source("../MainFunctions/vasicek_estimator.r")
source("../MainFunctions/al_omari_1_estimator.r")
source("../MainFunctions/al_omari_2_estimator.r")
source("../MainFunctions/bootstrap_estimator.r")

# genera muestras con repeticion:
bootstrap_estimator1 <- function(data, B, seed = NULL, verbose = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(data)
  t <- van_es_estimator(data)
  bx <- numeric(B)
  
  for (i in 1:B) {
    sampled_indices <- sample(1:n, n, replace = TRUE)
    sampled_data <- data[sampled_indices]
    bx[i] <- van_es_estimator(sampled_data)
    
    if (verbose) {
      cat("Bootstrap Sample", i, ": ", sampled_data, "\n")
    }
  }
 
  estimated_mean <- mean(bx)
  
  # cat("mean_entropies: ", estimated_mean, "\n")
  # cat("t: ", t, "\n")
  
  
  return(2*t - estimated_mean)
}

# seed1 <- set.seed(1234567890, kind = "Mersenne-Twister")
# # # Use
data <- c(3, 5, 1)
result <- bootstrap_estimator1(data, 5, seed = 123, verbose = TRUE)
 correa1 <- van_es_estimator(data)
#
 cat("entropy: ", result, "\n")
 cat("correa: ", correa1, "\n")
