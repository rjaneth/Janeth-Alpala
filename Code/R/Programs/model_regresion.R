# Clear the environment
rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggsci)
library(compiler)
library(minpack.lm)  # For nlsLM function


source("../../../Code/R/Programs/read_ENVI_images.R")

# Set seed for reproducibility
set.seed(1234567890, kind = "Mersenne-Twister")

# Load the SAR image data (assuming Z is your image matrix)
load("../Programs/Data/Phantom_4_z.Rdata")  # SAR image of 512x512 pixels
x <- Z  # Use x for consistency
rows <- nrow(x)
cols <- ncol(x)

# Parameters
L <- 5               # Number of looks
window_size <- 11    # Sliding window size
B <- 300            # Number of bootstrap samples within each window
sample.size <- window_size^2  # Number of pixels in each window
# Initialize matrix to store p-values of the nlsLM model
p_values_nlsLM <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterate over sliding windows
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Select local window
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
    window_vector <- as.vector(window_data)
    
    # Check if the window has enough variation
    if (length(unique(window_vector)) > 1) {
      # Initialize vectors to store bootstrapped means and CVs
      boot_means <- numeric(B)
      boot_cvs <- numeric(B)
      
      # Bootstrap resampling within the window
      for (b in 1:B) {
        resampled_data <- sample(window_vector, replace = TRUE)
        boot_means[b] <- mean(resampled_data)
        boot_cvs[b] <- sd(resampled_data) / mean(resampled_data)
      }
      
      # Prepare data for nlsLM
      data_boot <- data.frame(Mean = boot_means, CV = boot_cvs)
      
      # Fit the nlsLM model
      tryCatch({
        nls_model <- nlsLM(CV ~ sqrt(sample.size) * (1 - exp(-(beta0 + beta1 * Mean)) - beta2),
                           data = data_boot, start = list(beta0 = 0.1, beta1 = 0.01, beta2 = 0.01))
        
        # Extract the summary to get the p-value for beta1
        summary_nls <- summary(nls_model)
        p_value_beta1 <- coef(summary_nls)["beta1", "Pr(>|t|)"]  # p-value for beta1
        
        # Store the p-value
        p_values_nlsLM[i, j] <- p_value_beta1
      }, error = function(e) {
        # If the model fails to converge, store NA
        p_values_nlsLM[i, j] <- NA
      })
    } else {
      # If there's no variation in the window, store NA
      p_values_nlsLM[i, j] <- NA
    }
  }
  
  # Optional: Print progress
  if (i %% 10 == 0) {
    cat("Processed row", i, "out of", rows - window_size + 1, "\n")
  }
}

# Save the p-values matrix and the original image data
save(p_values_nlsLM, x, file = "./Data/p_values_nlsLM.Rdata")

# Load the p-values if necessary
 load("./Data/p_values_nlsLM.Rdata")

# Create a data frame for plotting
p_values_df <- melt(p_values_nlsLM, varnames = c("Row", "Col"), value.name = "p_value")

# Plot the p-values
ggplot(p_values_df, aes(x = Col, y = Row, fill = p_value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "plasma", na.value = "white", trans = "log", breaks = c(0.0001, 0.001, 0.01, 0.05, 1)) +
  labs(title = "P-values of beta1 Across the Image (nlsLM)",
       x = "Column",
       y = "Row",
       fill = "P-value") +
  theme_minimal() +
  coord_fixed()


heterogeneous_regions <- p_values_nlsLM < 0.05

heterogeneous_df <- melt(heterogeneous_regions, varnames = c("Row", "Col"), value.name = "Heterogeneous")

ggplot(heterogeneous_df, aes(x = Col, y = Row, fill = Heterogeneous)) +
  geom_raster() +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "blue", "NA" = "white")) +
  labs(title = "Identified Heterogeneous Regions",
       x = "Column",
       y = "Row",
       fill = "Heterogeneous") +
  theme_minimal() +
  coord_fixed()
