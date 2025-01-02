#445 boostrap

rm(list = ls())
source("../../../Code/R/MainFunctions/al_omari_1_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_al_omari_1_estimator.R")


source("../../../Code/R/MainFunctions/correa_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator_log_mean.R")


source("../../../Code/R/MainFunctions/ebrahimi_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_ebrahimi_estimator.R")

source("../../../Code/R/Programs/read_ENVI_images.R")
# small_image <- myread.ENVI(file='path/to/small_image_1.img', 
#                            headerfile='path/to/small_image_1.hdr')
# x <- myread.ENVI(file='../../../Data/SAR/Agua_envi/Intensity_HH.img', 
#                  headerfile='../../../Data/SAR/Agua_envi/Intensity_HH.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/forest_envi_100/Intensity_HH.img', 
#                  headerfile='../../../Data/SAR/forest_envi_100/Intensity_HH.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/crop2_envi/Intensity_HH.img', 
#                  headerfile='../../../Data/SAR/crop2_envi/Intensity_HH.hdr')
x <- myread.ENVI(file='../../../Data/SAR/urban_envi/Intensity_HH.img', 
                 headerfile='../../../Data/SAR/urban_envi/Intensity_HH.hdr')



# Función para calcular la estadística basada en entropía
calc_entropy_stat <- function(window_data, L, B) {
  return(bootstrap_al_omari_1_estimator(window_data, B) - 
           (log(mean(window_data)) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L))))
}

# # Leer la imagen pequeña
# small_image <- myread.ENVI(file='path/to/small_image.img', 
#                            headerfile='path/to/small_image.hdr')
set.seed(1234567890, kind = "Mersenne-Twister")
L <- 36
B <- 200
window_size <- 11
num_samples <- 100

# Preparar para almacenar resultados
entropy_stats <- numeric(num_samples)

# Tomar muestras con reposición
for (k in 1:num_samples) {
  # Seleccionar aleatoriamente una posición en la imagen pequeña
  i <- sample(1:(nrow(x) - window_size + 1), 1)
  j <- sample(1:(ncol(x) - window_size + 1), 1)
  
  # Seleccionar la ventana de 11x11 píxeles
  window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
  
  # Calcular la estadística de entropía para esta muestra
  entropy_stats[k] <- calc_entropy_stat(window_data, L, B)
}

# Calcular media y desviación estándar
mean_entropy_stat <- mean(entropy_stats, na.rm = TRUE)
sd_entropy_stat <- sd(entropy_stats, na.rm = TRUE)

# Guardar los resultados
save(entropy_stats, mean_entropy_stat, sd_entropy_stat, file = "./Data/results_entropy_urban_boot.Rdata")

# Repetir el proceso para las otras imágenes pequeñas (media y homogénea)
