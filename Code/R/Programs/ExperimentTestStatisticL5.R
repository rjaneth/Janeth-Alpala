library(ggplot2)
library(ggsci)


if(!require("rstudioapi")) install("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../MainFunctions/ebrahimi_estimator.R")
source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/correa_estimator.R")
source("../MainFunctions/bootstrap_correa_estimator.R")
source("../MainFunctions/bootstrap_correa_estimator_log_mean.R")


source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")






set.seed(1234567890, kind="Mersenne-Twister")

R <- 1000
mu <- 1
L <- 5
B <- 5

sample.size <- c(81)

TestStatistics <- NULL

for(s in sample.size){
  for(r in 1:R){
    
    z <- gamma_sar_sample(L, mu, s)
    
    TestStat <- bootstrap_correa_estimator(z, B) - log(mean(z))  + (-L + log(L) - lgamma(L) - (1 - L) * digamma(L)) 
    
    TestStatistics <- rbind(TestStatistics,
                            c(s, TestStat))
  }
  
  mean_entropy <- mean(TestStat)
  cat(" mean", mean_entropy, "\n" )
}

TestStatistics <- data.frame(TestStatistics)
names(TestStatistics) <- c("Sample Size", "Test Statistics")

#print(TestStatistics)


# ggplot(TestStatistics, aes(x=`Test Statistics`, col=`Sample Size`, group=`Sample Size`)) +
#    geom_density()

difference.betweeen.GammaSAR.GI0.OneLook <- function(alpha) {
  return(
    -lgamma(1-alpha) +
      (1-alpha) * digamma(1-alpha) -
      (1-alpha) * digamma(-alpha) +
      log(-1-alpha) +
      lgamma(-alpha) -1
  )
}


ggplot(TestStatistics, aes(x = `Test Statistics`, col = factor(`Sample Size`), linetype = factor(`Sample Size`))) +
  geom_line(stat = "density", linewidth = 1.0) +  
  scale_color_manual(
    values = pal_cosmic()(length(sample.size)),
    name = "Sample Size"
  ) +
  scale_linetype_manual(
    values = rep("solid", length(sample.size)),  
    name = "Sample Size"
  ) +
  # Expected values under H1
  geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.OneLook(-1.01), col="red") +
  annotate("text", x=difference.betweeen.GammaSAR.GI0.OneLook(-1.01), y=15, parse=TRUE, label="alpha==-1.01", hjust=-0.1, vjust=-1, col="red") +
  geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.OneLook(-1.1), col="red") +
  annotate("text", x=difference.betweeen.GammaSAR.GI0.OneLook(-1.1), y=15, parse=TRUE, label="alpha==-1.1", hjust=1.1, vjust=-1, col="red") +
  geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.OneLook(-1.5), col="red") +
  annotate("text", x=difference.betweeen.GammaSAR.GI0.OneLook(-1.5), y=15, parse=TRUE, label="alpha==-1.5", hjust=1.1, vjust=-1, col="red") +
  geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.OneLook(-300), col="red") +
  annotate("text", x=difference.betweeen.GammaSAR.GI0.OneLook(-300), y=15, parse=TRUE, label="alpha==-300", hjust=-.1, vjust=-1, col="red") +
  labs(x = "Test Statistics", y = "Density") +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    legend.position = "bottom",
    # legend.title = element_text(size = 10),
    # legend.text = element_text(size = 8),
    legend.key.size = unit(1, "lines")  
  )
