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

 x <- myread.ENVI(file='../../../Data/SAR/Flevoland_300/Intensity_VV.img', 
                  headerfile='../../../Data/SAR/Flevoland_300/Intensity_VV.hdr')

x <- myread.ENVI(file='../../../Data/SAR/Ottawa_512/Intensity_VV.img', 
                 headerfile='../../../Data/SAR/Ottawa_512/Intensity_VV.hdr')



L <- 5
#B <- 1




window_size <- 5


rows <- nrow(x)
cols <- ncol(x)




cv_values_cqv <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Seleccionar ventana local
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
    
    # Calcular los cuartiles
    q1 <- quantile(window_data, 0.25)
    q3 <- quantile(window_data, 0.75)
    
    # Calcular el CQV
    cv_values_cqv[i, j] <- (q3 - q1) / (q3 + q1)
  }
}

# Guardar el resultado en un archivo
save(cv_values_cqv, x, file = "./Data/results_data_Ottawa_cqv_512_5.Rdata")


