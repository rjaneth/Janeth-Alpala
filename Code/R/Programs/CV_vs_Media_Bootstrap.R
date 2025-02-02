#535
# Load necessary functions and libraries
rm(list = ls())
source("../MainFunctions/gi0_sample.R")
library(minpack.lm)

# Set parameters
#set.seed(1234567890, kind = "Mersenne-Twister")
sample.size <- 121
R <- 200
mu <- 1
L <- 5
alpha1 <- -2 # Very negative alpha

# Generate one sample from gi0 distribution
z <- gi0_sample(mu, alpha1, L, sample.size)

# Initialize vectors to store means and CVs
means <- numeric(R)
cvs <- numeric(R)

# Generate bootstrap samples and compute means and CVs with additional conditions
for (r in 1:R) {
  valid_sample <- FALSE
  while (!valid_sample) {
    # Generate bootstrap sample with replacement
    bootstrap_sample <- sample(z, size = sample.size, replace = TRUE)
    
    # Check if the sample has all unique values and finite results
    if (!all(bootstrap_sample == bootstrap_sample[1]) && 
        is.finite(mean(bootstrap_sample)) && 
        is.finite(sd(bootstrap_sample))) {
      
      # Calculate mean and CV
      means[r] <- mean(bootstrap_sample)
      cvs[r] <- sd(bootstrap_sample) / means[r]
      valid_sample <- TRUE
    }
  }
}

# Create data frame
data <- data.frame(Mean = means, CV = cvs)

# Plot the data to visualize the relationship
plot(data$Mean, data$CV, xlab = "Media", ylab = "CV", 
     main = paste("CV vs Media (Bootstrap), alpha =", alpha1), pch = 16, cex = 0.5)

# Adjust starting values based on data characteristics
# Since CV is nearly constant, we can simplify the model

# Simplified model: CV ~ constant
# Compute the average CV as the initial estimate
initial_cv <- mean(data$CV)

# Use nlsLM for robust fitting with adjusted starting values and bounds
nls_model_exp_simplified_lm <- nlsLM(
  CV ~ sqrt(sample.size) * (1 - exp(-(beta0 + beta1 * Mean)) - beta2),
  data = data,
  start = list(beta0 = -0.1, beta1 = 0.1, beta2 = 1),
  lower = c(-Inf, -Inf, -Inf),
  upper = c(Inf, Inf, Inf),
  control = nls.lm.control(maxiter = 1000)
)

summary(nls_model_exp_simplified_lm)

data$Fitted_CV_simplified_lm <- predict(nls_model_exp_simplified_lm)

# Plot the fitted curve
plot(data$Mean, data$CV, xlab = "Media", ylab = "CV", 
     main = paste("CV vs Media (Bootstrap), alpha =", alpha1), pch = 16, cex = 0.5)
lines(data$Mean[order(data$Mean)], 
      data$Fitted_CV_simplified_lm[order(data$Mean)], 
      col = "blue", lwd = 2)
