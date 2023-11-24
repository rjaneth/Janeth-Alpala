


 #source("../MainFunctions/van_es_estimator.r")
 
 bootstrap_van_es_estimator <- function(x, B){

   v.Bootstrap <- rep(0, B)

   for(b in 1:B) {
     same_values <- TRUE
     while (same_values) {
       bx <- sample(x, replace = TRUE)
       if (!all(bx == bx[1])) {
         same_values <- FALSE
         entropy_result <- van_es_estimator(bx)
         if (is.finite(entropy_result)) {
           v.Bootstrap[b] <- entropy_result
           #cat("Bootstrap Sample", b, ": ", bx, "\n")
           #cat("Entropy for Bootstrap Sample", b, ": ", v.Bootstrap[b], "\n")
         }
       }
     }
   }


   t <- van_es_estimator(x)
   estimated_mean <- mean(v.Bootstrap[!is.na(v.Bootstrap)])
   #cat("mean boot", estimated_mean, "\n")

   return(2*t - estimated_mean)

 }


# bootstrap_van_es_estimator <- function(data, B) {
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
#         entropy_result <- van_es_estimator(sampled_data)
#         if (is.finite(entropy_result)) {
#           bx[i] <- entropy_result
# 
#           #if (verbose) {
#             #cat("Bootstrap Sample", i, ": ", sampled_data, "\n")
#             #cat("Entropy for Bootstrap Sample", i, ": ", bx[i], "\n")
#           #}
#         }
#       }
#     }
#   }
# 
#   bx <- bx[!is.na(bx)]  # Elimina valores NA (si los hay)
#   estimated_mean <- mean(bx)
# 
#   return(2*t - estimated_mean)
# }

 # set.seed(1234567890, kind = "Mersenne-Twister")
 # data <-c(1.904591, 1.451404, 0.9608317, 1.229641, 0.429283, 0.6902405, 0.3525769, 0.4554712, 0.9444219)
 # result <- bootstrap_van_es_estimator(data, 5)
 # result1 <- van_es_estimator(data)
 # cat( "result :", result, "\n" )
 # cat( "result1 :", result1, "\n" )
