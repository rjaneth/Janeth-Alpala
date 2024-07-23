rm(list = ls())

library(ggplot2)
library(ggsci)
library(univariateML)
library(compiler)
library(MASS)
library(fitdistrplus)

#library(flexsurv)
#library(betaprime)
#library(actuar)
#library(MASS)
#library(fitdistrplus)
library(invgamma)
source("../MainFunctions/ebrahimi_estimator.R")
source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/correa_estimator.R")
source("../MainFunctions/bootstrap_correa_estimator.R")
source("../MainFunctions/bootstrap_correa_estimator_log_mean.R")
source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")

set.seed(1234567890, kind = "Mersenne-Twister")

sample.size <- c(49)
R <- 15000
mu <- 1
L <- 5
B <- 5
alpha1 <- -3

TestStatistics <- numeric(0)  

for (s in sample.size) {
  for (r in 1:R) {
    #z <- gi0_sample(mu, alpha1, L, s)
    z <- gamma_sar_sample(L, mu, s)
    TestStat <- sd(z) / mean(z)
    TestStatistics <- c(TestStatistics, TestStat)  
  }
}

iqr <- IQR(TestStatistics)


n <- length(TestStatistics)


bins_fd <- (2 * iqr) / (n^(1/3))


num_bins_fd <- round(diff(range(TestStatistics)) / bins_fd)


print(num_bins_fd)
# iqr <- IQR(TestStatistics)
# 
# # Número de datos
# n <- length(TestStatistics)
# 
# # Aplicar la regla de Freedman-Diaconis
# bins_fd <- (2 * iqr) / (n^(1/3))
# 
# 
# num_bins_fd <- round(diff(range(TestStatistics) / bins_fd)
# 
# 
# print(num_bins_fd)

fln <- fitdist(TestStatistics, "lnorm")
fg <- fitdist(TestStatistics, "gamma")
fn <- fitdist(TestStatistics, "norm")  
fw <- fitdist(TestStatistics, "weibull")

# 
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot.legend <- c( "Lognormal", "Gamma", "Normal", "Weibull") 
num_bins <-86

denscomp(list( fln, fg, fn, fw), n = num_bins, legendtext = plot.legend, fitcol = c("darkblue", "#637029","indianred3","darkgrey"), fitlwd=c(3,3,3,3), fitlty = 1,xlab = "CV ", xlegend = "topright", ylegend = NULL)
legend("topright", legend = plot.legend, col = c("darkblue", "#637029","indianred3","darkgrey"), lwd = 3, lty = 1, box.lwd = 1,box.col = "white",bg = "white")
# denscomp(list( fln, fg, fn,fw), legendtext = plot.legend, fitcol = c("darkblue", "#637029","indianred3","darkgrey"), fitlwd=c(3,3,3,3), fitlty = 1,xlab = "CV ")
qqcomp(list(fln, fg, fn,fw), legendtext = plot.legend,fitcol = c("darkblue", "#637029","indianred3","darkgrey"), fitlwd=c(3,3,3,3),xlegend = "topleft", ylegend = NULL)
cdfcomp(list(fln, fg, fn,fw), legendtext = plot.legend, fitcol = c("darkblue", "#637029","indianred3","darkgrey"), fitlwd=c(3,3,3,3), fitlty = 1,xlab = "CV")
ppcomp(list(fln, fg, fn,fw), legendtext = plot.legend,fitcol = c("darkblue", "#637029","indianred3","darkgrey"), fitlwd=c(3,3,3,3),xlegend = "topleft", ylegend = NULL)
