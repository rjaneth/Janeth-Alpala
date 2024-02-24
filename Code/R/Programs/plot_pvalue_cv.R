
load("./Data/results_data_Flevoland_cv_300_5.Rdata")
alpha <- 5.452
beta <- 2.061


pgamma_inv <- function(x, alpha, beta) {
  1 - pgamma(1/x, shape = alpha, rate = beta)
}


p_values_matrix <- matrix(NA, nrow = nrow(cv_values), ncol = ncol(cv_values))

for (i in 1:nrow(cv_values)) {
  for (j in 1:ncol(cv_values)) {
    cv_value <- cv_values[i, j]
    
    
    p_values_matrix[i, j] <- pgamma_inv(cv_value, alpha, beta)
  }
}


save(p_values_matrix, file = "./Data/results_pvalue_cv_Flevoland_cv_300_5.Rdata")


source("../imagematrix.R")
plot(imagematrix(p_values_matrix > 0.05))
