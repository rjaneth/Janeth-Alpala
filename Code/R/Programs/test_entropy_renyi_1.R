# 492
# Start the timer
#rm(list = ls())
start_time <- Sys.time()

#Estimators
source("../../../Code/R/MainFunctions/al_omari_1_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_al_omari_1_estimator.R")


source("../../../Code/R/MainFunctions/renyi_entropy_estimator_v1.R")
source("../../../Code/R/MainFunctions/bootstrap_renyi_entropy_estimator_v1.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator_log_mean.R")


source("../../../Code/R/MainFunctions/ebrahimi_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_ebrahimi_estimator.R")

source("../../../Code/R/Programs/read_ENVI_images.R")




# x <- myread.ENVI(file='../../../Data/SAR/panama512/Intensity_VV.img', 
#                   headerfile='../../../Data/SAR/panama512/Intensity_VV.hdr')

#  x <- myread.ENVI(file='../../../Data/SAR/mexico_512/Intensity_VV.img', 
#                   headerfile='../../../Data/SAR/mexico_512/Intensity_VV.hdr')
# # 

#x <- myread.ENVI(file='../../../Data/SAR/Rotterdam_1024/Intensity_HH.img', 
#                 headerfile='../../../Data/SAR/Rotterdam_1024/Intensity_HH.hdr')
#
x <- myread.ENVI(file='../../../Data/SAR/chicago_1024/Intensity_HH.img', 
                 headerfile='../../../Data/SAR/chicago_1024/Intensity_HH.hdr')
#

rows <- nrow(x)
cols <- ncol(x)



L <- 36 
B <- 100 
lambda<- 0.98
window_size <- 7 


difference_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterate over sliding windows
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Select local window
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]#change for x if the data comes SAR images
       difference_values[i, j] <-bootstrap_renyi_entropy_estimator_v1(window_data, B, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
                                                                                        (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(window_data))-log(L))
    #difference_values[i, j] <-renyi_entropy_estimator_v1(window_data, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
    #                                                                               (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(window_data))-log(L))
    
    # difference_values[i, j] <-bootstrap_renyi_entropy_estimator_v1(window_data, B, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
    #                                                                                             (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(window_data))-log(L))
    #                                                         ((lambda * L - lambda + 1) * log(lambda)) -   lambda * lgamma(L))+log(mean(window_data)))
  }
}



save(difference_values, x, file = "./Data/results_Illinois_2_renyi_1024_w7_b100_L36_098.Rdata")
# Stop the timer
end_time <- Sys.time()


# Calculate the elapsed time
execution_time <- end_time - start_time
execution_time

