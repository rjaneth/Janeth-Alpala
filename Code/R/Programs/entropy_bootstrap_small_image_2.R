#447 file
# Función para calcular la estadística basada en entropía
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


x <- myread.ENVI(file='../../../Data/SAR/urban_envi/Intensity_HH.img', 
                 headerfile='../../../Data/SAR/urban_envi/Intensity_HH.hdr')

calc_entropy_stat <- function(window_data, L, B) {
  return(bootstrap_al_omari_1_estimator(window_data, B) - 
           (log(mean(window_data)) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L))))
}

# # Leer una de las imágenes pequeñas
# small_image <- myread.ENVI(file='path/to/small_image.img', 
#                            headerfile='path/to/small_image.hdr')

L <- 36
B <- 200
window_size <- 11
num_samples <- 1000
sample_size <- window_size * window_size

# Preparar para almacenar resultados
entropy_stats <- numeric(num_samples)

# Tomar muestras con reposición
for (k in 1:num_samples) {
  # Seleccionar una muestra aleatoria de tamaño 121
  i <- sample(1:(nrow(x) - window_size + 1), 1)
  j <- sample(1:(ncol(x) - window_size + 1), 1)
  window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
  
  # Calcular la estadística de entropía para esta muestra
  entropy_stats[k] <- calc_entropy_stat(window_data, L, B)
}

# Calcular media y desviación estándar
mean_entropy_stat <- mean(entropy_stats, na.rm = TRUE)
sd_entropy_stat <- sd(entropy_stats, na.rm = TRUE)

# Guardar los resultados
save(entropy_stats, mean_entropy_stat, sd_entropy_stat, file = "./Data/results_entropy_small_image.Rdata")
