# Start the timer
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
# Load simulated image EX. Phantom
#load("../Programs/Data/Phantom_4_regions.Rdata")
#save(Z, file="./Data/Phantom_4_regions.Rdata")
load("../Programs/Data/Phantom_4_z.Rdata")

# Assign the dimensions of the loaded matrix to the variables rows and cols
# change for x if the data comes SAR images

rows <- nrow(Z)
cols <- ncol(Z)

# rows <- nrow(x)
# cols <- ncol(x)



L <- 5 # Number of looks
B <- 200 # Replications bootstrap
lambda<- 1.3
window_size <- 7 # sliding window size


difference_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterate over sliding windows
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Select local window
    window_data <- Z[i:(i + window_size - 1), j:(j + window_size - 1)]#change for x if the data comes SAR images
    difference_values[i, j] <-renyi_entropy_estimator_v1(window_data, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
                                                                                                (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(window_data))-log(L))
   #   difference_values[i, j] <-bootstrap_renyi_entropy_estimator_v1(window_data, B, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
   #                                                                                     (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(window_data))-log(L))
   # # #apply the entropy test using one of the bootstrap estimators above
    #difference_values[i, j] <- bootstrap_al_omari_1_estimator(window_data,B) - (log(mean(window_data)) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L)))
    #difference_values[i, j] <-bootstrap_renyi_entropy_estimator(window_data, B, lambda) - (-log(L)+ (1 / (1 - lambda)) * (lgamma(lambda * L - lambda + 1) - 
                           # bootstrap_renyi_entropy_estimator_v1(z, B, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
                                                                        #(lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(z))-log(L))                                                            ((lambda * L - lambda + 1) * log(lambda)) -   lambda * lgamma(L))+log(mean(window_data)))
  }
}


# Save the results
#save(difference_values, x, file = "./Data/results_Chicago_1024_9_AO_200b_36L.Rdata")
#save(difference_values, Z, file = "./Data/results_Phantom_4_renyi_B100_w7_075_L5.Rdata")#w9
##save(difference_values, Z, file = "./Data/results_Phantom_4_renyi_B100.Rdata")
#save(difference_values, Z, file = "./Data/results_Phantom_4_renyi_B200_w9_08_L5.Rdata")
#save(difference_values, Z, file = "./Data/results_Phantom_4_renyi_B100_w7_09_L5.Rdata")
#save(difference_values, Z, file = "./Data/results_Phantom_4_AO_B100_w7_L5.Rdata")#ok
#save(difference_values, Z, file = "./Data/results_Phantom_4_AO_B100_w9_L5.Rdata")
#save(difference_values, Z, file = "./Data/results_Phantom_4_renyi_B200_w7_085_L5.Rdata") 
#save(difference_values, Z, file = "./Data/results_Phantom_4_renyi_w7_09_L5.Rdata")#OK
save(difference_values, Z, file = "./Data/results_Phantom_4_renyi_w7_1_3_L5.Rdata")#1_6, 2_0
# Stop the timer
end_time <- Sys.time()


# Calculate the elapsed time
execution_time <- end_time - start_time
execution_time

