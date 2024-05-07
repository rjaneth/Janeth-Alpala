rm(list = ls())

# if(!require("rstudioapi")) install("rstudioapi")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../../../Code/R/MainFunctions/bootstrap_correa_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator_log_mean.R")
source("../../../Code/R/Programs/read_ENVI_images.R")
source("../../../Code/R/MainFunctions/correa_estimator.R")

# x <- myread.ENVI(file='../../../Data/SAR/chicago_1024/Intensity_HH.img', 
#                  headerfile='../../../Data/SAR/chicago_1024/Intensity_HH.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/mexico_512/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/mexico_512/Intensity_VV.hdr')

x <- myread.ENVI(file='../../../Data/SAR/Lake_512/Intensity_VV.img', 
                 headerfile='../../../Data/SAR/Lake_512/Intensity_VV.hdr')



rows <- nrow(x)
cols <- ncol(x)

window_size <- 9


MnADmedian <- function(mat){
  mdn <- median(mat)
  MnAD <- mean(abs(mat - mdn))
  MnAD/mdn
}

cd_values_mnad <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Seleccionar ventana local
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
    
    # Calcular el CD basado en MnAD para la ventana local (matriz)
    cd_values_mnad[i, j] <- MnADmedian(window_data)
  }
}

# Guardar el resultado en un archivo
save( cd_values_mnad, x, file = "./Data/MnAD_results_data_lake_512.Rdata")
#save( cd_values_mnad, x, file = "./Data/MnAD_results_data_mexico_512.Rdata")
#save( cd_values_mnad, x, file = "./Data/MnAD_results_data_Illinois_crops_1024.Rdata")
#save(cd_values_mnad, x, file = "./Data/results_data_Ottawa_512_mnad_7.Rdata")
