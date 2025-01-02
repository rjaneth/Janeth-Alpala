rm(list = ls())
#449
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

# x <- myread.ENVI(file='../../../Data/SAR/Lake_512/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/Lake_512/Intensity_VV.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/Rotterdam_1024/Intensity_HH.img', 
#                  headerfile='../../../Data/SAR/Rotterdam_1024/Intensity_HH.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/urban_envi/Intensity_HH.img', 
#                  headerfile='../../../Data/SAR/urban_envi/Intensity_HH.hdr')

small_image <- myread.ENVI(file='../../../Data/SAR/forest_envi_100/Intensity_HH.img', 
                 headerfile='../../../Data/SAR/forest_envi_100/Intensity_HH.hdr')

# Función para calcular el estadístico MnAD basado en la mediana
calc_mnad_stat <- function(window_data) {
  mdn <- median(window_data)
  MnAD <- mean(abs(window_data - mdn))
  return(MnAD / mdn)
}

# Leer la imagen pequeña
# small_image <- myread.ENVI(file='path/to/small_image.img', 
#                            headerfile='path/to/small_image.hdr')

window_size <- 11
num_samples <- 100

# Preparar para almacenar resultados
mnad_stats <- numeric(num_samples)

# Tomar muestras con reposición
for (k in 1:num_samples) {
  # Seleccionar aleatoriamente una posición en la imagen pequeña
  i <- sample(1:(nrow(small_image) - window_size + 1), 1)
  j <- sample(1:(ncol(small_image) - window_size + 1), 1)
  
  # Seleccionar la ventana de 9x9 píxeles
  window_data <- small_image[i:(i + window_size - 1), j:(j + window_size - 1)]
  
  # Calcular el estadístico MnAD basado en la mediana para esta muestra
  mnad_stats[k] <- calc_mnad_stat(window_data)
}

# Calcular media y desviación estándar
mean_mnad_stat <- mean(mnad_stats, na.rm = TRUE)
sd_mnad_stat <- sd(mnad_stats, na.rm = TRUE)

# Guardar los resultados
#save(mnad_stats, mean_mnad_stat, sd_mnad_stat, file = "./Data/results_mnad_small_image.Rdata")
