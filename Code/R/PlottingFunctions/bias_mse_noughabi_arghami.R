rm(list = ls())
library(ggplot2)
library(gtools)
library(stats4)
library(rmutil)
#library(gsl)
library(invgamma)

library(ggthemes)
library(grid)
library(gridExtra)
library(gtable)
#install.packages("cowplot")
library(cowplot)
library(ggpubr)
theme_set(theme_bw()  +
            theme(text=element_text(family="serif"),
                  legend.position = "top")# Gridtop , right , bottom , or left#, panel.grid = element_blank()
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

source("../MainFunctions/bootstrap_van_es_estimator.r")
source("../MainFunctions/bootstrap_correa_estimator.r")
source("../MainFunctions/bootstrap_ebrahimi_estimator.r")
source("../MainFunctions/bootstrap_noughabi_arghami_estimator.r")
source("../MainFunctions/bootstrap_vasicek_estimator.r")
source("../MainFunctions/bootstrap_al_omari_1_estimator.r")
source("../MainFunctions/bootstrap_al_omari_2_estimator.r")

set.seed(1234567890, kind="Mersenne-Twister")

sample_sizes <- c(9, 25, 49, 81, 121)
R <- 10
mu <- 1
L <- 1

# 
results <- data.frame()

for (n in sample_sizes) {
  true_entropy <- entropy_gamma_sar(L, mu)
  bias_est <- rep(0, R)
  mse_est <- rep(0, R)
  bias_bootstrap <- rep(0, R)
  mse_bootstrap <- rep(0, R)
  
  for (r in 1:R) {
    # Samples
    data <- gamma_sar_sample(L, mu, n)
    
    # Entropy estimations
    est <- noughabi_arghami_estimator(data)
    bootstrap_est <- bootstrap_noughabi_arghami_estimator(data, B = 10)
    
    #  bias and MSE
    bias_est[r] <- mean(est) - true_entropy
    mse_est[r] <- (est - true_entropy)^2
    bias_bootstrap[r] <- mean(bootstrap_est) - true_entropy
    mse_bootstrap[r] <- (bootstrap_est - true_entropy)^2
  }
  
  # 
  result_row <- data.frame(
    Sample_Size = n,
    Estimator = c("Noughabi-Arghami", "Bootstrap Noughabi-Arghami"),
    Bias = c(mean(bias_est), mean(bias_bootstrap)),
    MSE = c(mean(mse_est), mean(mse_bootstrap))
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
