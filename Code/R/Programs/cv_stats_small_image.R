#448

# Código para leer una imagen pequeña y aplicar el test estadístico
rm(list = ls())
# Leer una de las imágenes pequeñas
# Supongamos que las imágenes pequeñas están guardadas como 'small_image_1', 'small_image_2', 'small_image_3'
source("../../../Code/R/MainFunctions/al_omari_1_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_al_omari_1_estimator.R")


source("../../../Code/R/MainFunctions/correa_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator_log_mean.R")


source("../../../Code/R/MainFunctions/ebrahimi_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_ebrahimi_estimator.R")

source("../../../Code/R/Programs/read_ENVI_images.R")

# Función para calcular la estadística basada en el coeficiente de variación (CV)
calc_cv_stat <- function(window_data) {
  mean_values <- mean(window_data)
  sd_values <- sd(window_data)
  return(sd_values / mean_values)
}
 # x <- myread.ENVI(file='../../../Data/SAR/Agua_envi/Intensity_HH.img', 
 #                  headerfile='../../../Data/SAR/Agua_envi/Intensity_HH.hdr')
  # x <- myread.ENVI(file='../../../Data/SAR/forest_envi_100/Intensity_HH.img', 
  #                 headerfile='../../../Data/SAR/forest_envi_100/Intensity_HH.hdr')
x <- myread.ENVI(file='../../../Data/SAR/urban_envi/Intensity_HH.img', 
                 headerfile='../../../Data/SAR/urban_envi/Intensity_HH.hdr')


window_size <- 11
num_samples <- 100

# Preparar para almacenar resultados
cv_stats <- numeric(num_samples)

# Tomar muestras con reposición
for (k in 1:num_samples) {
  # Seleccionar aleatoriamente una posición en la imagen pequeña
  i <- sample(1:(nrow(x) - window_size + 1), 1)
  j <- sample(1:(ncol(x) - window_size + 1), 1)
  
  # Seleccionar la ventana de 9x9 píxeles
  window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
  
  # Calcular la estadística de CV para esta muestra
  cv_stats[k] <- calc_cv_stat(window_data)
}

# Calcular media y desviación estándar
mean_cv_stat <- mean(cv_stats, na.rm = TRUE)
sd_cv_stat <- sd(cv_stats, na.rm = TRUE)

# Guardar los resultados
save(cv_stats, mean_cv_stat, sd_cv_stat, file = "./Data/results_cv_small_image.Rdata")
