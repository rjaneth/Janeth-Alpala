#rm(list = ls())

# if(!require("rstudioapi")) install("rstudioapi")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../../../Code/R/MainFunctions/bootstrap_correa_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator_log_mean.R")
source("../../../Code/R/Programs/read_ENVI_images.R")
source("../../../Code/R/MainFunctions/correa_estimator.R")

#x <- myread.ENVI(file='../../../Data/SAR/Ottawa/Intensity_VV.img', headerfile='../../../Data/SAR/Ottawa/Intensity_VV.hdr')
# x_Flevoland2 <- myread.ENVI(file='../../../Data/SAR/Flevoland_100/Intensity_VV.img', 
#                             headerfile='../../../Data/SAR/Flevoland_100/Intensity_VV.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/Flevoland_300/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/Flevoland_300/Intensity_VV.hdr')

x <- myread.ENVI(file='../../../Data/SAR/Ottawa_512/Intensity_VV.img', 
                 headerfile='../../../Data/SAR/Ottawa_512/Intensity_VV.hdr')



L <- 5
#B <- 1




window_size <- 5


rows <- nrow(x)
cols <- ncol(x)


cv_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)


#test_difference_vector <- numeric()

# 
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Seleccionar ventana local
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
    
     
    x_values <- window_data
    
    
    # 
    mean_values <- mean(window_data)
    sd_values <- sd(window_data)
    
    # 
    cv_values[i, j] <- sd_values / mean_values
    

   
  }
}


#save(test_difference_vector, file = "./Data/results_data_Flevoland_100_5.Rdata")
save( cv_values, x, file = "./Data/results_data_Ottawa_cv_512_5.Rdata")



# print(mean_values)
# print(entropy_values)
# print(true_entropy_values)
# print(difference_values)
# print(test_difference_vector)

# # Load data from R.data file
# load("results_data.Rdata")
# 
# # Display information
# cat("Mean Values:\n")
# print(mean_values)
# 
# cat("Entropy Values:\n")
# print(entropy_values)
# 
# cat("True Entropy Values:\n")
# print(true_entropy_values)
# 
# cat("Difference Values:\n")
# print(difference_values)
# 
# cat("Test Difference Vector:\n")
# print(test_difference_vector)
