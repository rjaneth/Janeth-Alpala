library(ggplot2)
library(reshape2)
#library(plotly)
library(knitr)
library(pandoc)
library(gridExtra)
library(gtools)
library(stats4)
library(rmutil)
library(scales)
library(tidyr)
library(rmutil)
library(invgamma)
library(tidyverse)
#library(RColorBrewer)
library(ggsci)
#library(wesanderson)
library(ggpubr)
library(patchwork)
library(dplyr)
#options(kableExtra.latex.load_packages = FALSE)
library(devtools)
#devtools::install_github("haozhu233/kableExtra")
library(kableExtra)
library(ggthemes)
library(latex2exp)

theme_set(theme_minimal() +
            theme(text=element_text(family="serif"),
                  legend.position = "bottom")#  top , right , bottom , or left#, panel.grid = element_blank()
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



# 
sample_sizes <- c(9, 25, 49, 81)
R <- 10
mu <- 10
L <- 5
alpha <- -100

#  Gamma SAR
true_entropy <- entropy_gamma_sar(L, mu)
#

generate_samples <- function(sample_size, replication, mu, alpha, L) {
  samples <- vector("list", replication)
  for (r in 1:replication) {
    samples[[r]] <- gi0_sample(mu, alpha, L, sample_size)#
  }
  return(samples)
}

 
calculate_entropy_and_test_table <- function(sample_sizes, R, mu, alpha, L) {
  estimators <- list(
    "Correa Bootstrap" = bootstrap_correa_estimator
  )
  
  output <- data.frame(
    SampleSize = integer(0),
    Estimator = character(0),
    MeanEntropy = numeric(0),
    ZStatistic = numeric(0),
    PValue = numeric(0)
  )
  
  for (ssize in sample_sizes) {
    samples <- generate_samples(ssize, R, mu, alpha, L)
    
    for (estimator_name in names(estimators)) {
      estimator <- estimators[[estimator_name]]
      v.entropy <- numeric(R)
      
      for (r in 1:R) {
        sample <- samples[[r]]
        
        if (grepl("Bootstrap", estimator_name)) {
          v.entropy[r] <- estimator(sample, B = 30)
        } else {
          v.entropy[r] <- estimator(sample)
        }
      }
      
      mean_entropy <- mean(v.entropy)
      
      # 
      z_statistic <- sqrt(R) * (mean_entropy - true_entropy) / sd(v.entropy)
      p_value <- 2 * (1 - pnorm(abs(z_statistic)))
      
      output <- rbind(
        output,
        data.frame(
          SampleSize = ssize,
          Estimator = estimator_name,
          MeanEntropy = round(mean_entropy, 5),
          ZStatistic = round(z_statistic, 5),
          PValue = round(p_value, 5)
        )
      )
    }
  }
  
  return(output)  
}

result_table <- calculate_entropy_and_test_table(sample_sizes, R, mu, alpha, L)
print(result_table)


 