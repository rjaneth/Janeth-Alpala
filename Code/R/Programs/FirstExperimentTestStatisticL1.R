library(ggplot2)
library(ggsci)


if(!require("rstudioapi")) install("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../MainFunctions/ebrahimi_estimator.R")
source("../MainFunctions/gamma_sar_sample.R")




set.seed(1234567890, kind="Mersenne-Twister")

R <- 500
mu <- 1
sample.size <- c(49, 81, 121)

TestStatistics <- NULL

for(s in sample.size){
  for(r in 1:R){
    
    z <- gamma_sar_sample(1, mu, s)
    TestStat <- ebrahimi_estimator(z) - (1 + log(mean(z)))
    #TestStat <-sqrt(R) *(ebrahimi_estimator(z) - (1 + log(mean(z))))/ sd(z)
    #cat(" test2", TestStat, "\n")
    TestStatistics <- rbind(TestStatistics,
                            c(s, TestStat))
  }
  mean_test <- mean(TestStat)
  cat(" mean", mean_test, "\n")
}

TestStatistics <- data.frame(TestStatistics)
names(TestStatistics) <- c("Sample Size", "Test Statistics")

# ggplot(TestStatistics, aes(x=`Test Statistics`, col=`Sample Size`, group=`Sample Size`)) +
#    geom_density()

difference.betweeen.GammaSAR.GI0.OneLook <- function(alpha) {
  return(
    -log(gamma(1-alpha)) +
      (1-alpha) * digamma(1-alpha) -
      (1-alpha) * digamma(-alpha) +
      log(-1-alpha) +
      log(gamma(-alpha)) -1
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
  #Expected values under H1
  geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.OneLook(-1.01), col="red") +
  annotate("text", x=difference.betweeen.GammaSAR.GI0.OneLook(-1.01), y=15, parse=TRUE, label="alpha==-1.01", hjust=-0.1, vjust=-1, col="red") +
  geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.OneLook(-1.1), col="red") +
  annotate("text", x=difference.betweeen.GammaSAR.GI0.OneLook(-1.1), y=15, parse=TRUE, label="alpha==-1.1", hjust=1.1, vjust=-1, col="red") +
  geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.OneLook(-1.5), col="red") +
  annotate("text", x=difference.betweeen.GammaSAR.GI0.OneLook(-1.5), y=15, parse=TRUE, label="alpha==-1.5", hjust=1.1, vjust=-1, col="red") +
  geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.OneLook(-3), col="red") +
  annotate("text", x=difference.betweeen.GammaSAR.GI0.OneLook(-3), y=15, parse=TRUE, label="alpha==-3", hjust=-.1, vjust=-1, col="red") +
  labs(x = "Test Statistics", y = "Density") +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    legend.position = "bottom",
    # legend.title = element_text(size = 10),
    # legend.text = element_text(size = 8),
    legend.key.size = unit(1, "lines")  
  )
