# Clear the environment
rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(mgcv)        # For GAMs
library(reshape2)
library(gridExtra)
library(ggsci)
library(compiler)

# Source your custom functions
source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/entropy_gamma_sar.R")
source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")
source("../../../Code/R/Programs/read_ENVI_images.R")

# Set seed for reproducibility
set.seed(1234567890, kind = "Mersenne-Twister")
# Load the SAR image data (assuming Z is your image matrix)
load("../Programs/Data/Phantom_4_z.Rdata") # SAR image of 512x512 pixels
x <- Z  # Use x for consistency
rows <- nrow(x)
cols <- ncol(x)

# Parameters
L <- 5               # Number of looks
window_size <- 11    # Sliding window size
B <- 100             # Number of bootstrap samples within each window
sample.size <- window_size^2  # Number of pixels in each window
# Initialize matrix to store p-values of the smooth term
p_values_smooth <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

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
      
      # Prepare data for GAM
      data_boot <- data.frame(Mean = boot_means, CV = boot_cvs)
      
      # Fit the GAM
      tryCatch({
        gam_model <- gam(CV ~ s(Mean), data = data_boot, family = gaussian())
        
        # Extract p-value of the smooth term
        summary_gam <- summary(gam_model)
        p_value_smooth <- summary_gam$s.pv[1]  # p-value for the smooth term
        
        # Store the p-value
        p_values_smooth[i, j] <- p_value_smooth
      }, error = function(e) {
        # If the model fails to converge, store NA
        p_values_smooth[i, j] <- NA
      })
    } else {
      # If there's no variation in the window, store NA
      p_values_smooth[i, j] <- NA
    }
  }
  # Optional: Print progress
  if (i %% 10 == 0) {
    cat("Processed row", i, "out of", rows - window_size + 1, "\n")
  }
}

# Save the p-values matrix and the original image data
save(p_values_smooth, x, file = "./Data/p_values_smooth_gam.Rdata")

# Load the p-values if necessary
# load("./Data/p_values_smooth_gam.Rdata")

# Create a data frame for plotting
p_values_df <- melt(p_values_smooth, varnames = c("Row", "Col"), value.name = "p_value")

# Plot the p-values
ggplot(p_values_df, aes(x = Col, y = Row, fill = p_value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "plasma", na.value = "white", trans = "log", breaks = c(0.0001, 0.001, 0.01, 0.05, 1)) +
  labs(title = "P-values of Smooth Term Across the Image (GAM)",
       x = "Column",
       y = "Row",
       fill = "P-value") +
  theme_minimal() +
  coord_fixed()

# Threshold the p-values to identify regions
# For example, regions where p-value < 0.05 are considered heterogeneous
heterogeneous_regions <- p_values_smooth < 0.05

# Plot the heterogeneous regions
heterogeneous_df <- melt(heterogeneous_regions, varnames = c("Row", "Col"), value.name = "Heterogeneous")

ggplot(heterogeneous_df, aes(x = Col, y = Row, fill = Heterogeneous)) +
  geom_raster() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue", "NA" = "white")) +
  labs(title = "Identified Heterogeneous Regions (GAM)",
       x = "Column",
       y = "Row",
       fill = "Heterogeneous") +
  theme_minimal() +
  coord_fixed()
