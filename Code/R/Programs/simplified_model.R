#493
# Clear the environment
rm(list = ls())

# Load necessary libraries
library(ggplot2)
library(minpack.lm)  # For robust nonlinear regression
library(univariateML)
library(compiler)
library(MASS)
library(fitdistrplus)
library(invgamma)
library(gridExtra)

# Load any additional libraries you need
library(ggsci)

# Source your custom functions
source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/entropy_gamma_sar.R")
source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")
source("../../../Code/R/Programs/read_ENVI_images.R")


# Set seed for reproducibility
set.seed(1234567890, kind = "Mersenne-Twister")

# Load and read data from SAR images
# x <- myread.ENVI(
#   file = '../../../Data/SAR/Rotterdam_1024/Intensity_HH.img',
#   headerfile = '../../../Data/SAR/Rotterdam_1024/Intensity_HH.hdr'
# )

load("../Programs/Data/Phantom_4_z.Rdata") # imagen SAR de 512x512 pixeles

rows <- nrow(Z)
cols <- ncol(Z)
#rows <- nrow(x)
#cols <- ncol(x)

# Parameters
L <- 5  # Number of looks
window_size <- 81  # Sliding window size
B <- 100  # Number of bootstrap samples within each window
sample.size <- window_size^2  # Number of pixels in each window

# Initialize matrix to store p-values of beta1
p_values_beta1 <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterate over sliding windows
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Select local window
    window_data <- Z[i:(i + window_size - 1), j:(j + window_size - 1)]
    window_vector <- as.vector(window_data)
    
    # Initialize vectors to store bootstrapped means and CVs
    boot_means <- numeric(B)
    boot_cvs <- numeric(B)
    
    # Check if the window has enough variation
    if (length(unique(window_vector)) > 1) {
      # Bootstrap resampling within the window
      for (b in 1:B) {
        resampled_data <- sample(window_vector, replace = TRUE)
        boot_means[b] <- mean(resampled_data)
        boot_cvs[b] <- sd(resampled_data) / mean(resampled_data)
      }
      
      # Prepare data for regression
      data_boot <- data.frame(Mean = boot_means, CV = boot_cvs)
      
      # Starting values for beta0 and beta1
      start_values <- list(beta0 = 0.1, beta1 = 0.1)
      
      # Fit the regression model
      # We may need to fix parameters or use a simplified model due to limited data 
      #model <- nls(CV ~ sqrt(sample.size) * (1 - exp(-(beta0 + beta1 * Mean)))
      tryCatch({
        model <- nlsLM(
          CV ~ sqrt(sample.size) * (1 - exp(-(beta0 + beta1 * Mean))),
          data = data_boot,
          start = start_values,
          control = nls.lm.control(maxiter = 1000)
        )
        
        # Extract p-value of beta1
        summary_model <- summary(model)
        p_value_beta1 <- coef(summary_model)[2, "Pr(>|t|)"]  # beta1 is the second parameter
        
        # Store the p-value
        p_values_beta1[i, j] <- p_value_beta1
      }, error = function(e) {
        # If the model fails to converge, store NA
        p_values_beta1[i, j] <- NA
      })
    } else {
      # If there's no variation in the window, store NA
      p_values_beta1[i, j] <- NA
    }
  }
  # Optional: Print progress
  if (i %% 10 == 0) {
    cat("Processed row", i, "out of", rows - window_size + 1, "\n")
  }
}

# Save the p-values matrix and the original image data
save(p_values_beta1, Z, file = "./Data/p_values_beta1_v1_11_200.Rdata")

# Optionally, you can visualize the p-values matrix
# For example, plot regions where p-value is less than 0.05
library(reshape2)

# Create a data frame for plotting
#load("./Data/p_values_beta1_v1.Rdata")
load("./Data/p_values_beta1_v1_11.Rdata")

p_values_df <- melt(p_values_beta1, varnames = c("Row", "Col"), value.name = "p_value")

# Plot the p-values
ggplot(p_values_df, aes(x = Col, y = Row, fill = p_value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "plasma", na.value = "white") +
  labs(title = "P-values of Beta1 Across the Image",
       x = "Column",
       y = "Row",
       fill = "P-value") +
  theme_minimal() +
  coord_fixed()

# You can also threshold the p-values to identify regions
# For example, regions where p-value < 0.05 are considered heterogeneous
heterogeneous_regions <- p_values_beta1 < 0.05

# Plot the heterogeneous regions
heterogeneous_df <- melt(heterogeneous_regions, varnames = c("Row", "Col"), value.name = "Heterogeneous")

ggplot(heterogeneous_df, aes(x = Col, y = Row, fill = Heterogeneous)) +
  geom_raster() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue", "NA" = "white")) +
  labs(title = "Identified Heterogeneous Regions",
       x = "Column",
       y = "Row",
       fill = "Heterogeneous") +
  theme_minimal() +
  coord_fixed()
