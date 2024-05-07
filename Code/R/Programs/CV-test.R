rm(list = ls())

# if(!require("rstudioapi")) install("rstudioapi")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# source("../../../Code/R/MainFunctions/bootstrap_correa_estimator.R")
# source("../../../Code/R/MainFunctions/bootstrap_correa_estimator_log_mean.R")
 source("../../../Code/R/Programs/read_ENVI_images.R")
# source("../../../Code/R/MainFunctions/correa_estimator.R")


# x <- myread.ENVI(file='../../../Data/SAR/chicago_1024/Intensity_HH.img', 
#                   headerfile='../../../Data/SAR/chicago_1024/Intensity_HH.hdr')

 # x <- myread.ENVI(file='../../../Data/SAR/mexico_512/Intensity_VV.img', 
 #                  headerfile='../../../Data/SAR/mexico_512/Intensity_VV.hdr')
 
  # x <- myread.ENVI(file='../../../Data/SAR/Lake_512/Intensity_VV.img', 
  #                  headerfile='../../../Data/SAR/Lake_512/Intensity_VV.hdr')
  load("../Programs/Data/Phantom_4_z.Rdata")#Phantom_4_z
  # Asignar las dimensiones de la matriz cargada a las variables rows y cols
  rows <- nrow(Z)
  cols <- ncol(Z)
  
  window_size <- 7
#L <- 5
#B <- 100


window_size <- 9


rows <- nrow(x)
cols <- ncol(x)


cv_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)



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

save( cv_values, x, file = "./Data/CV_results_data_lake_512.Rdata")
#save( cv_values, x, file = "./Data/CV_results_data_Mexixo_512.Rdata")
#save( cv_values, x, file = "./Data/CV_results_data_Illinois_crops_1024.Rdata")
#save(cv_values, Z,  file = "./Data/results_data_simulated_Phantom_7.Rdata")


