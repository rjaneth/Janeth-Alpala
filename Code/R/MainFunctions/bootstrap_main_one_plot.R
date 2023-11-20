rm(list = ls())
library(ggplot2)
library(tidyr)
library(gtools)
library(stats4)
library(rmutil)

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
source("../MainFunctions/bootstrap_van_es_estimator.r")
source("../MainFunctions/bootstrap_correa_estimator.r")
source("../MainFunctions/bootstrap_ebrahimi_estimator.r")
source("../MainFunctions/bootstrap_noughabi_arghami_estimator.r")
source("../MainFunctions/bootstrap_vasicek_estimator.r")
source("../MainFunctions/bootstrap_al_omari_1_estimator.r")
source("../MainFunctions/bootstrap_al_omari_2_estimator.r")

set.seed(1234567890, kind = "Mersenne-Twister")
sample_sizes <- c(9, 25, 49)
#sample_sizes <- c(9, 25)
R <-5
mu <- 1
L <- 1


output <- NULL


for(ssize in sample_sizes) {
  v.bootstrap.entropy <- NULL
  v.nonparametric.entropy <- NULL
  for(r in 1:R) {
    sample <- gamma_sar_sample(L, mu, ssize)
    cat("muestra: ", sample, "\n")
    v.bootstrap.entropy[r] <- bootstrap_van_es_estimator(sample, B=3)
    cat("boot: ", v.bootstrap.entropy, "\n")
    #v.bootstrap.entropy[r] <- est_bootstrap
    v.nonparametric.entropy[r] <- van_es_estimator(sample)
    cat("van es: ", v.nonparametric.entropy, "\n")
  }
  output <- rbind(output, data.frame(SampleSize = ssize,  Nonparametric = mean(v.nonparametric.entropy), Bootstrap = mean(v.bootstrap.entropy)))
}
print(output)
#cat("entropy: ", output, "\n")




output <- gather(output, key = "Estimator", value = "Mean", -SampleSize)

ggplot(output, aes(x = SampleSize, y = Mean, color = Estimator)) +
 geom_hline(yintercept = entropy_gamma_sar(1, 1)) +
  geom_line() +
  geom_point() +
  xlab("Sample Size") +
  ylab("Mean Entropy") +
  #scale_color_manual(values = c("Bootstrap" = "blue", "Nonparametric" = "red")) +
  #scale_color_manual(values = scales::rainbow_pal()(2)) + 
  labs(color = "Estimator") +  
  theme_minimal()


# output <- data.frame(output)
# names(output) <- c("Mean_estimator", "Mean_bootstrap")
# output$SampleSize <- sample_sizes
# #plot
# ggplot(output, aes(x=SampleSize, y=Mean)) +
#   geom_hline(yintercept = entropy_gamma_sar(1, 1)) +
#   geom_line() +
#   geom_point() +
#   #geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.1) +
#   xlab("Sample Size") +
#   #scale_x_log10(breaks=output$SampleSize) +
#   ylab("Mean Non-parametric Entropy")