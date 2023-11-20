rm(list = ls())
library(ggplot2)
library(gtools)
#library(stats4)
#library(rmutil)
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


set.seed(1234567890, kind = "Mersenne-Twister")



# Set sample sizes, the number of repetitions, and parameters
sample_sizes <- c(9, 25, 49, 81, 121, 225)
R <- 1000
mu <- 1
L <- 1


# Define the true entropy value for the gamma SAR distribution
true_entropy <- entropy_gamma_sar(L, mu)  # Replace with the true entropy value
# Initialize the output data frame
output <- NULL

# Perform non-parametric entropy estimation for different sample sizes
for (ssize in sample_sizes) {
  v.nonparametric.entropy <- NULL
  for (r in 1:R) {
    sample <- gamma_sar_sample(L, mu, ssize)
    v.nonparametric.entropy[r] <- van_es_estimator(sample)
  }
  
  # Calculate the bias by subtracting the true entropy value
  bias <- mean(v.nonparametric.entropy) - true_entropy
  
  # Store the bias in the output data frame
  output <- rbind(output, c(bias, ssize))
}

# Create a data frame with results
output <- data.frame(output)
names(output) <- c("Bias")
output$SampleSize <- sample_sizes

# Plot the bias as a function of the sample size
ggplot(output, aes(x=SampleSize, y=Bias)) +
  geom_line() +
  geom_point() +
  xlab("Sample Size") +
  ylab("Bias of Non-parametric Entropy Estimator")
