#537
# Load necessary functions and libraries
rm(list = ls())
source("../MainFunctions/gi0_sample.R")
library(minpack.lm)

# Set parameters
#set.seed(1234567890, kind = "Mersenne-Twister")
sample.size <- 200
R <- 3000
mu <- 1
L <- 5
alpha1 <- -1.5 # Very negative alpha

# Generate one sample from gi0 distribution
z <- gi0_sample(mu, alpha1, L, sample.size)

# Initialize vectors to store means and CVs
means <- numeric(R)
cvs <- numeric(R)

# Generate bootstrap samples with variable size and additional conditions
for (r in 1:R) {
  valid_sample <- FALSE
  while (!valid_sample) {
    # Generate a variable bootstrap sample size within 80% to 100% of original size
    bootstrap.size <- sample(seq(round(sample.size * 0.95), sample.size), 1)
    bootstrap_sample <- sample(z, size = bootstrap.size, replace = TRUE)
    
    # Calculate mean and CV
    current_mean <- mean(bootstrap_sample)
    current_cv <- sd(bootstrap_sample) / current_mean
    
    # Check conditions: variability in mean and CV, finite results, and sufficient diversity
    if (!all(bootstrap_sample == bootstrap_sample[1]) && 
        is.finite(current_mean) && is.finite(current_cv) &&
        (r == 1 || abs(current_mean - means[r-1]) > 0.05) && 
        (r == 1 || abs(current_cv - cvs[r-1]) > 0.05)) {
      
      means[r] <- current_mean
      cvs[r] <- current_cv
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
  start = list(beta0 = -0.1, beta1 = 0.1, beta2 = 0),
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
