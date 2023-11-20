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

source("../MainFunctions/van_es_estimator.r")
source("../MainFunctions/correa_estimator.r")
source("../MainFunctions/ebrahimi_estimator.r")
source("../MainFunctions/noughabi_arghami_estimator.r")
source("../MainFunctions/vasicek_estimator.r")
source("../MainFunctions/al_omari_1_estimator.r")
source("../MainFunctions/al_omari_2_estimator.r")

source("../MainFunctions/bootstrap_estimator.r")
source("../MainFunctions/bootstrap_van_es_estimator.r")
source("../MainFunctions/bootstrap_correa_estimator.r")
source("../MainFunctions/bootstrap_ebrahimi_estimator.r")
source("../MainFunctions/bootstrap_noughabi_arghami_estimator.r")
source("../MainFunctions/bootstrap_vasicek_estimator.r")
source("../MainFunctions/bootstrap_al_omari_1_estimator.r")
source("../MainFunctions/bootstrap_al_omari_2_estimator.r")

set.seed(1234567890, kind = "Mersenne-Twister")
sample_sizes <- c(9, 25,49, 81, 121)
R <- 100
mu <- 1
L <- 1

true_entropy <- entropy_gamma_sar(L, mu)
#output <- NULL
results <- data.frame()
for(ssize in sample_sizes) {
  v.bootstrap.entropy <- NULL
  v.nonparametric.entropy <- NULL
  
  for(r in 1:R) {
    sample <- gamma_sar_sample(L, mu, ssize)
   
    v.bootstrap.entropy[r] <- bootstrap_van_es_estimator(sample, B = 100) 
    
    v.nonparametric.entropy[r] <- van_es_estimator(sample)
    
  }
  
  mse_nonparametric <- mean((v.nonparametric.entropy - true_entropy)^2)
  
  mse_bootstrap <- mean((v.bootstrap.entropy - true_entropy)^2)
  
  
  bias_nonparametric <- mean(v.nonparametric.entropy) - true_entropy
  
  bias_bootstrap <- mean(v.bootstrap.entropy) - true_entropy
  
  
  result_row <- data.frame(
    Sample_Size = ssize,
    Estimator = c("Van Es", "Bootstrap Van Es"),
    Bias = c((bias_nonparametric), (bias_bootstrap)),
    MSE = c(mse_nonparametric, mse_bootstrap)
  )
  results <- rbind(results, result_row)
}





print(results)


bias_plot <- ggplot(results, aes(x = Sample_Size, y = Bias, color = Estimator)) +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_point() +
  labs(title = "",
       y = "Bias", x = "Sample size") 


# Plot MSE
mse_plot <- ggplot(results, aes(x = Sample_Size, y = MSE, color = Estimator)) +
  geom_hline(yintercept = 0) +
  geom_line() +
  geom_point() +
  labs(title = "",
       y = "MSE", x = "Sample size") 

ggarrange(bias_plot, mse_plot, ncol=2, nrow=1, common.legend = TRUE, legend="top")


