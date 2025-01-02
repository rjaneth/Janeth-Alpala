#Codigo que usa la entropia de shannon en el estimador y test
# Start the timer
start_time <- Sys.time()

#Estimators
source("../../../Code/R/MainFunctions/al_omari_1_estimator.R") #estimador no parametrico
source("../../../Code/R/MainFunctions/bootstrap_al_omari_1_estimator.R")# estimador no parametrico mejorado, que vamos a usar

source("../../../Code/R/Programs/read_ENVI_images.R")


#Load and read data from SAR images

# x <- myread.ENVI(file='../../../Data/SAR/panama512/Intensity_VV.img', 
#                   headerfile='../../../Data/SAR/panama512/Intensity_VV.hdr')

load("../Programs/Data/Phantom_4_z.Rdata") #imagen simulada

# Assign the dimensions of the loaded matrix to the variables rows and cols
# change for x if the data comes SAR images

rows <- nrow(Z)
cols <- ncol(Z)

# rows <- nrow(x)
# cols <- ncol(x)



L <- 5 # Number of looks
B <- 100 # Replications bootstrap
window_size <- 7 # sliding window size


difference_values_shannon <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterate over sliding windows
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Select local window
    window_data <- Z[i:(i + window_size - 1), j:(j + window_size - 1)]#change for x if the data comes SAR images
    difference_values_shannon[i, j] <-renyi_entropy_estimator_v1(window_data, lambda) -((lambda * lgamma(L) - lgamma(lambda * (L - 1) + 1) +
                                                                                   (lambda * (L - 1) + 1) * log(lambda)) / (lambda - 1) + log(mean(window_data))-log(L))
 
  }
}


# Save the results

#save(difference_values_shannon, Z, file = "./Data/results_Phantom_4_AO_B100_w7_L5.Rdata")

# Stop the timer
end_time <- Sys.time()


# Calculate the elapsed time
execution_time <- end_time - start_time
execution_time

