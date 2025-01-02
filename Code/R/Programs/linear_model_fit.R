# Load necessary libraries
rm(list = ls())
library(ggplot2)

# Source your custom functions
source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/entropy_gamma_sar.R")
source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")

# Set seed for reproducibility
set.seed(1234567890, kind = "Mersenne-Twister")

# Simulation parameters
sample.size <- 144
R <- 10000
mu <- 1
L <- 5
alpha1 <- -1.5 # You can change this value for different scenarios

# Initialize vectors to store means and CVs
means <- numeric(R)
cvs <- numeric(R)

# Simulate data
for (r in 1:R) {
  z <- gi0_sample(mu, alpha1, L, sample.size)
  means[r] <- mean(z)
  cvs[r] <- sd(z) / mean(z)
}

# Prepare data frame
data <- data.frame(Mean = means, CV = cvs)

# Calculate U and Z
sqrt_n <- sqrt(sample.size)
data$U <- data$CV / sqrt_n

# Remove observations where U >= 1 to avoid taking log of zero or negative numbers
data <- subset(data, U < 1)

# Calculate Z
data$Z <- -log(1 - data$U)

# Fit the linear model
model_linear <- lm(Z ~ Mean, data = data)

# Summary of the model
summary(model_linear)

# Extract coefficients
coefficients <- coef(model_linear)

# Predict Z for plotting
mean_seq <- seq(min(data$Mean), max(data$Mean), length.out = 1000)
predicted_Z <- predict(model_linear, newdata = data.frame(Mean = mean_seq))

# Calculate predicted CV from predicted Z
predicted_U <- 1 - exp(-predicted_Z)
predicted_CV <- predicted_U * sqrt_n

# Create a data frame for the predicted values
prediction_df <- data.frame(Mean = mean_seq, CV = predicted_CV)

# Plot the data and the fitted model
ggplot(data, aes(x = Mean, y = CV)) +
  geom_point(alpha = 0.3, color = "red") +
  geom_line(data = prediction_df, aes(x = Mean, y = CV), color = "blue", size = 1) +
  labs(title = "Scatter Plot of CV versus Mean with Linear Model Fit",
       x = "Mean",
       y = "Coefficient of Variation (CV)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90")
  )
