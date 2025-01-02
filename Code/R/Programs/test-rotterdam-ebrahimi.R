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


#Load and read data from SAR images

# x <- myread.ENVI(file='../../../Data/SAR/panama512/Intensity_VV.img', 
#                   headerfile='../../../Data/SAR/panama512/Intensity_VV.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/chicago_1024/Intensity_HH.img', 
#                  headerfile='../../../Data/SAR/chicago_1024/Intensity_HH.hdr')
# 

 y <- myread.ENVI(file='../../../Data/SAR/Rotterdam_1024/Intensity_HH.img',
                  headerfile='../../../Data/SAR/Rotterdam_1024/Intensity_HH.hdr')

#y <- myread.ENVI(file='../../../Data/SAR/envi_Rotterdam_2_1024/Intensity_HH.img',
#                 headerfile='../../../Data/SAR/envi_Rotterdam_2_1024/Intensity_HH.hdr')#new

# x <- myread.ENVI(file='../../../Data/SAR/envi_panama_512/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/envi_panama_512/Intensity_VV.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/chicago_1024/Intensity_HH.img', 
#                  headerfile='../../../Data/SAR/chicago_1024/Intensity_HH.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/SanFrancisco_650/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/SanFrancisco_650/Intensity_VV.hdr')

# Load simulated image EX. Phantom
#load("../Programs/Data/Phantom_4_regions.Rdata")
#save(Z, file="./Data/Phantom_4_regions.Rdata")
#load("../Programs/Data/Phantom_4_z.Rdata")

# Assign the dimensions of the loaded matrix to the variables rows and cols
# change for x if the data comes SAR images

#rows <- nrow(Z)
#cols <- ncol(Z)

rows <- nrow(y)
cols <- ncol(y)



#illinois 1024
# L <- 36 # Number of looks
# B <- 100 # Replications bootstrap
# lambda<- 0.85
# window_size <- 9 # sliding window size
# L <- 36 # Number of looks
# B <- 200 # Replications bootstrap
# lambda<- 0.75
# window_size <- 11 # sliding window size

L <- 1 # Number of looks
#B <- 200 # Replications bootstrap
lambda<- 0.99
window_size <- 7 # sliding window size

difference_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterate over sliding windows
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Select local window
    window_data <- y[i:(i + window_size - 1), j:(j + window_size - 1)]#change for x if the data comes SAR images
    # difference_values[i, j] <-renyi_entropy_estimator_v1(window_data, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
    #                                                                                (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(window_data))-log(L))
    #apply the entropy test using one of the bootstrap estimators above
    difference_values[i, j] <- ebrahimi_estimator(window_data) - (log(mean(window_data)) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L)))
    #difference_values[i, j] <-bootstrap_renyi_entropy_estimator(window_data, B, lambda) - (-log(L)+ (1 / (1 - lambda)) * (lgamma(lambda * L - lambda + 1) - 
    # bootstrap_renyi_entropy_estimator_v1(z, B, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
    #(lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(z))-log(L))                                                            ((lambda * L - lambda + 1) * log(lambda)) -   lambda * lgamma(L))+log(mean(window_data)))
  }
}


# Save the results
#save(difference_values, x, file = "./Data/results_Chicago_1024_9_AO_200b_36L.Rdata")
#save(difference_values, Z, file = "./Data/results_Phantom_4_renyi_B100.Rdata")
#save(difference_values, x, file = "./Data/results_rotterdam_renyi_w7_b200_0_65.Rdata")
#save(difference_values, x, file = "./Data/results_panama_512_renyi_w7_b100_075_L5_w9.Rdata")
#save(difference_values, x, file = "./Data/results_illinois_1024_renyi_b100_085_L36_w9.Rdata")
#save(difference_values, x, file = "./Data/results_sanFrancisco_650_renyi_b200_075_L5_w11.Rdata")
#save(difference_values, x, file = "./Data/results_illinois_1024_renyi_b200_075_L36_w11.Rdata")
#save(difference_values, x, file = "./Data/results_illinois_1024_renyi_b200_09_L36_w7.Rdata")
save(difference_values, y, file = "./Data/results_rotterdam_origin_Ebrahimi_w7.Rdata")
# Stop the timer
end_time <- Sys.time()


# Calculate the elapsed time
execution_time <- end_time - start_time
execution_time

