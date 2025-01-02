# Load necessary libraries
rm(list = ls())
library(ggplot2)
library(ggsci)
library(MASS)  # For Gamma GLM
library(gridExtra)

# Source your custom functions
source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/entropy_gamma_sar.R")
source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")

# Set seed for reproducibility
set.seed(1234567890, kind = "Mersenne-Twister")

# Simulation parameters
sample.size <- 49
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

# Fit the Gamma GLM with log link
model_glm <- glm(CV ~ Mean, family = Gamma(link = "log"), data = data)

# Summary of the model
summary(model_glm)

# Extract coefficients
coefficients <- coef(model_glm)

# Create a sequence of Mean values for plotting the fitted curve
mean_seq <- seq(min(data$Mean), max(data$Mean), length.out = 1000)
predicted_cv <- predict(model_glm, newdata = data.frame(Mean = mean_seq), type = "response")

ggplot(data, aes(x = Mean, y = CV)) +
  geom_point(alpha = 0.3, color = "red") +
  geom_line(
    data = data.frame(Mean = mean_seq, CV = predicted_cv),
    aes(x = Mean, y = CV),
    color = "blue",
    size = 1
  ) +
  labs(
    title = "Scatter Plot of CV versus Mean with Gamma GLM Fit",
    x = "Mean",
    y = "Coefficient of Variation (CV)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90")
  )
