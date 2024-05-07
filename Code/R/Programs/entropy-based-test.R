# Start the timer
start_time <- Sys.time()

#Estimators
source("../../../Code/R/MainFunctions/al_omari_1_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_al_omari_1_estimator.R")


source("../../../Code/R/MainFunctions/correa_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator.R")
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
window_size <- 7 # sliding window size


difference_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterate over sliding windows
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Select local window
    window_data <- Z[i:(i + window_size - 1), j:(j + window_size - 1)]#change for x if the data comes SAR images
   #apply the entropy test using one of the bootstrap estimators above
    difference_values[i, j] <- bootstrap_al_omari_1_estimator(window_data,B) - (log(mean(window_data)) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L)))
    
  }
}

# Save the results
#save(difference_values, x, file = "./Data/results_Chicago_1024_9_AO_200b_36L.Rdata")
save(difference_values, Z, file = "./Data/results_Phantom_4_z1_200b.Rdata")

# Stop the timer
end_time <- Sys.time()


# Calculate the elapsed time
execution_time <- end_time - start_time
execution_time

