
rm(list = ls())
if(!require("rstudioapi")) install("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(e1071)
library(nortest)

load("./Data/results_data_Flevoland_100.Rdata")



df <- data.frame(TestDifference = test_difference_vector)


ggplot(df, aes(x = TestDifference)) +
  geom_density(color = "blue", fill = "lightblue", alpha = 0.7) +
  labs(x = "Test Difference", y = "Density") +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    legend.position = "bottom"
  )



# 
mean_val <- mean(test_difference_vector)
sd_val <- sd(test_difference_vector)
var_val <- var(test_difference_vector)
skewness_val <- skewness(test_difference_vector)
kurtosis_val <- kurtosis(test_difference_vector)
ad_p_value <- ad.test(test_difference_vector)$p.value
cv_val <- sd_val / mean_val

#
cat("Mean:", mean_val, "\n")
cat("Standard Deviation:", sd_val, "\n")
cat("Variance:", var_val, "\n")
cat("Coefficient of Variation:", cv_val, "\n")
cat("Skewness:", skewness_val, "\n")
cat("Kurtosis:", kurtosis_val, "\n")
cat("Anderson-Darling p-value:", ad_p_value, "\n")

# 
summary_stats <- data.frame(
  Mean = mean_val,
  SD = sd_val,
  Variance = var_val,
  CV = cv_val,
  Skewness = skewness_val,
  Kurtosis = kurtosis_val,
  adpvalue = ad_p_value
)

qq_plot <- qqnorm(test_difference_vector, plot.it = FALSE)
qq_data <- data.frame(Theoretical = qq_plot$x, Sample = qq_plot$y)


ggplot(qq_data, aes(x = Theoretical, y = Sample)) +
  geom_point() +
  geom_abline(intercept = mean(test_difference_vector), slope = sd(test_difference_vector), col = "red", lty = 2) +
  labs(title = "QQ Plot", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()


save(summary_stats, qq_data, file = "./Data/results.Rdata")


calculate_p_value <- function(test_difference, mu, sigma) {
  epsilon <- (test_difference - mu) / sigma
  p_value <- 2 * (1 - pnorm(abs(epsilon)))
  
  return(p_value)
}

# abs_test_difference_vector <- abs(test_difference_vector)
# p_values <- sapply(abs_test_difference_vector, function(value) {
#   calculate_p_value(value, mean_val, sd_val)
# })


Calcular p-values para cada valor en test_difference_vector
p_values <- sapply(test_difference_vector, function(value) {
  calculate_p_value(value, mean_val, sd_val)
})


p_values_df <- data.frame(
  TestDifference = abs_test_difference_vector,
  PValue = p_values
)





# ggplot(p_values_df, aes(x = PValue, y = TestDifference)) +
#   geom_point(color = "blue", size = 2, alpha = 0.7) +
#   labs(title = " ", x = "P-value", y = "Test Difference") +
#   theme_minimal()