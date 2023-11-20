
# This code demonstrates non-parametric entropy estimation for Gamma SAR data.
# It uses non-parametric estimators to calculate the sample mean and standard deviation 
# of the entropy for various sample sizes.
# The code also includes error bars based on the standard deviation for visualization.


rm(list = ls())
library(ggplot2)

source("../MainFunctions/gamma_sar_sample.r")
source("../MainFunctions/van_es_estimator.r")
source("../MainFunctions/entropy_gamma_sar.r")
set.seed(1234567890, kind="Mersenne-Twister")

# Define sample sizes, the number of repetitions, and parameters
sample_sizes <- c(9, 25, 49, 81, 121, 225)
R <- 100
mu <- 1
L <- 1
true_entropy <- entropy_gamma_sar(2, 2)
#cat( "True entropy:", tru2_entropy)
# Initialize the output data frame
output <- NULL

# Perform non-parametric entropy estimation for different sample sizes
for(ssize in sample_sizes) {
  v.nonparametric.entropy <- NULL
  for(r in 1:R) {
    sample <- gamma_sar_sample(L, mu, ssize)
    v.nonparametric.entropy[r] <- van_es_estimator(sample)
  }
  output <- rbind(output, c(mean(v.nonparametric.entropy), sd(v.nonparametric.entropy)))
  cat( "Sample Size:", ssize, "entropy:", mean(v.nonparametric.entropy),"True entropy:", true_entropy, "\n")
}

# Create a data frame with results
output <- data.frame(output)
names(output) <- c("Mean", "Std")
output$SampleSize <- sample_sizes

#cat( "True entropy:", tru2_entropy)
#Plot
ggplot(output, aes(x=SampleSize, y=Mean)) +
  geom_hline(yintercept = entropy_gamma_sar(1, 1)) +
  geom_line() +
  geom_point() +
  #geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.1) +
  xlab("Sample Size") +
  #scale_x_log10(breaks=output$SampleSize) +
  ylab("Mean Non-parametric Entropy")

# Save the plot as a PDF file
#ggsave("myplot.pdf")


