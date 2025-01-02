#490
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

 x <- myread.ENVI(file='../../../Data/SAR/Rotterdam_1024/Intensity_HH.img',
                  headerfile='../../../Data/SAR/Rotterdam_1024/Intensity_HH.hdr')
 
 # x <- myread.ENVI(file='../../../Data/SAR/envi_panama_512/Intensity_VV.img', 
 #                  headerfile='../../../Data/SAR/envi_panama_512/Intensity_VV.hdr')
 
 # x <- myread.ENVI(file='../../../Data/SAR/chicago_1024/Intensity_HH.img', 
 #                  headerfile='../../../Data/SAR/chicago_1024/Intensity_HH.hdr')
 
 # x <- myread.ENVI(file='../../../Data/SAR/SanFrancisco_650/Intensity_VV.img',
 #                  headerfile='../../../Data/SAR/SanFrancisco_650/Intensity_VV.hdr')
 # x <- myread.ENVI(file='../../../Data/SAR/Amazonas_700/Intensity_HH.img',
 #                  headerfile='../../../Data/SAR/Amazonas_700/Intensity_HH.hdr')

 # x <- myread.ENVI(file='../../../Data/SAR/Lake_512/Intensity_VV.img', 
 #                  headerfile='../../../Data/SAR/Lake_512/Intensity_VV.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/envi_sevilla_2500/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/envi_sevilla_2500/Intensity_VV.hdr')
# x <- myread.ENVI(file='../../../Data/SAR/sinaloa_512/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/sinaloa_512/Intensity_VV.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/Ottawa_512_N/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/Ottawa_512_N/Intensity_VV.hdr')

#rows <- nrow(Z)
#cols <- ncol(Z)

 rows <- nrow(x)
 cols <- ncol(x)

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
 B <- 50 # Replications bootstrap
 lambda<- 0.7
 window_size <- 7 # sliding window size

difference_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterate over sliding windows
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Select local window
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]#change for x if the data comes SAR images
    difference_values[i, j] <-bootstrap_renyi_entropy_estimator_v1(window_data, B, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
                                                                                                (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(window_data))-log(L))
    #apply the entropy test using one of the bootstrap estimators above
    #difference_values[i, j] <- bootstrap_al_omari_1_estimator(window_data,B) - (log(mean(window_data)) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L)))
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
#save(difference_values, x, file = "./Data/results_illinois1_512_renyi_b200_08_L36_w9.Rdata")
#save(difference_values, x, file = "./Data/results_spain_2500_renyi_b100_08_L5_w9.Rdata")
#save(difference_values, x, file = "./Data/results_sinaloa_512_renyi_b100_085_L18_w9.Rdata")
#save(difference_values, x, file = "./Data/results_sinaloa_512_renyi_b100_09_L18_w7.Rdata")
#save(difference_values, x, file = "./Data/results_Ottawa_512_N_renyi_b100_085_L5_w7.Rdata")
#save(difference_values, x, file = "./Data/results_SanFrancisco_650_renyi_b100_085_L5_w7.Rdata")#no funciona
#save(difference_values, x, file = "./Data/results_SanFrancisco_650_renyi_b1_099_L1_w9.Rdata")
#save(difference_values, x, file = "./Data/results_SanFranciscoN_400_renyi_b200_08_L5_w9.Rdata")# no funciona
#save(difference_values, x, file = "./Data/results_Amazonas_700_renyi_b100_085_L18_w7.Rdata")# no funciona
save(difference_values, x, file = "./Data/results_Rotterdam_1024_renyi_b50_07_L1_w7.Rdata")# no funciona

# Stop the timer
end_time <- Sys.time()


# Calculate the elapsed time
execution_time <- end_time - start_time
execution_time

