# Non-Parametric Shannon Entropy Estimator
# Description: This function calculates Shannon entropy using a non-parametric estimator,
# known as the van Es estimator. 
# It is based on spacings between data points 
#
# Input:
#   - data: sample or dataset 
#
# Output:
#   - Estimated Shannon entropy 

van_es_estimator <- function(data) {
  n <- length(data)
  m <- round(sqrt(n) + 0.5)  # m-spacing
  data_sorted <- sort(data)
  sum_term1 <- 0
  sum_term2 <- 0
  
  for (i in 1:(n - m)) {
    sum_term1 <- sum_term1 + log(((n + 1) / m) * (data_sorted[i + m] - data_sorted[i]))
  }
  
  for (k in m:n) {
    sum_term2 <- sum_term2 + 1 / k
  }
  
  sum_term3 <- log(m / (n + 1))
  
  return((sum_term1 / (n - m)) + sum_term2+ sum_term3)
  
}

# # Example usage:
# data <- c(1.2, 1.5,3.6, 2.0, 2.1, 2.5, 2.6, 3.0, 3.1, 3.5 )
# result <- van_es_estimator(data)
# cat("Avan es :", result, "\n")