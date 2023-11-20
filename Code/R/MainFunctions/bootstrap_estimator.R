


# source("../MainFunctions/van_es_estimator.r")
# source("../MainFunctions/correa_estimator.r")
# source("../MainFunctions/ebrahimi_estimator.r")
# source("../MainFunctions/noughabi_arghami_estimator.r")
# source("../MainFunctions/vasicek_estimator.r")
# source("../MainFunctions/al_omari_1_estimator.r")
# source("../MainFunctions/al_omari_2_estimator.r")
# #source("../MainFunctions/bootstrap_estimator.r")

bootstrap_estimator <- function(data, B, verbose = TRUE) {
  n <- length(data)
  t <- van_es_estimator(data)
  bx <- numeric(B)
  
  for (i in 1:B) {
    same_values <- TRUE
    while (same_values) {
      sampled_indices <- sample(1:n, n, replace = TRUE)
      sampled_data <- data[sampled_indices]
      if (!all(sampled_data == sampled_data[1])) {
        same_values <- FALSE
        entropy_result <- van_es_estimator(sampled_data)
        if (is.finite(entropy_result)) {
          bx[i] <- entropy_result
          
          if (verbose) {
            #cat("Bootstrap Sample", i, ": ", sampled_data, "\n")
            #cat("Entropy for Bootstrap Sample", i, ": ", bx[i], "\n")
          }
        }
      }
    }
  }
  
  bx <- bx[!is.na(bx)]  # Elimina valores NA (si los hay)
  estimated_mean <- mean(bx)
  
  return(2*t - estimated_mean)
}

# bootstrap_estimator <- function(data, B, verbose = TRUE) {
#   #if (!is.null(seed)) set.seed(seed)
#   n <- length(data)
#   t <- van_es_estimator(data)
#   bx <- numeric(B)
#   
#   for (i in 1:B) {
#     same_values <- TRUE
#     while (same_values) {
#       sampled_indices <- sample(1:n, n, replace = TRUE)
#       sampled_data <- data[sampled_indices]
#       if (!all(sampled_data == sampled_data[1])) {
#         same_values <- FALSE
#         bx[i] <- van_es_estimator(sampled_data)
#         
#         if (verbose) {
#           cat("Bootstrap Sample", i, ": ", sampled_data, "\n")
#           cat("Entropy for Bootstrap Sample", i, ": ", bx[i], "\n")
#         }
#       }
#     }
#   }
#   #hist(bx, main = "Histograma de bx", xlab = "Valores de bx", ylab = "Frecuencia", col = "lightblue")
#   
#   estimated_mean <- mean(bx)
#   
#   # cat("mean_entropies: ", estimated_mean, "\n")
#   # cat("t: ", t, "\n")
#   
#   return(2*t - estimated_mean)
# }

# # seed1 <- set.seed(1234567890, kind = "Mersenne-Twister")
# # # # Use
# data <- c(3, 5, 1, 8, 2, 7, 4, 9, 6, 3)
# result <- bootstrap_Estimator(data, 1000, seed = 123, verbose = TRUE)
# correa1 <- correa_estimator(data)
# 
# cat("entropy: ", result, "\n")
# cat("correa: ", correa1, "\n")