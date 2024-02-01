rm(list = ls())

if(!require("rstudioapi")) install("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("../../../Code/R/MainFunctions/bootstrap_correa_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator_log_mean.R")
source("../../../Code/R/Programs/read_ENVI_images.R")
source("../../../Code/R/MainFunctions/correa_estimator.R")

#x <- myread.ENVI(file='../../../Data/SAR/Ottawa/Intensity_VV.img', headerfile='../../../Data/SAR/Ottawa/Intensity_VV.hdr')
x_Flevoland2 <- myread.ENVI(file='../../../Data/SAR/Flevoland_100/Intensity_VV.img', 
                            headerfile='../../../Data/SAR/Flevoland_100/Intensity_VV.hdr')

L <- 5
B <- 1

#difference_total <- bootstrap_correa_estimator_log_mean(x_Flevoland2, B)-  (L - log(L) + lgamma(L) + (1 - L) * digamma(L)) 

# Tamaño de la ventana
window_size <- 15


rows <- nrow(x_Flevoland2)
cols <- ncol(x_Flevoland2)

# 
mean_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)
entropy_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)
true_entropy_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)
difference_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)


test_difference_vector <- numeric()

# Iterar sobre ventanas deslizantes
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Seleccionar ventana local
    window_data <- x_Flevoland2[i:(i + window_size - 1), j:(j + window_size - 1)]

    # Almacenar valores de x para cada ventana 
    x_values <- window_data

    # Calcular la media
    mean_values[i, j] <- mean(window_data)

    # Calcular entropía con el estimador no paramétrico
    entropy_values[i, j] <- correa_estimator(window_data)

    # Calcular entropía verdadera
    true_entropy_values[i, j] <- (log(mean(window_data)) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L)))

    # Calcular la diferencia entre entropías
    
    difference_values[i, j] <- bootstrap_correa_estimator_log_mean(window_data, B)-  (L - log(L) + lgamma(L) + (1 - L) * digamma(L)) 
      #correa_estimator(window_data)- (log(mean(window_data)) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L)))
    test_difference_vector <- c(test_difference_vector, difference_values[i, j])
  }
}

# Save  vector to R.data file
save(test_difference_vector, file = "./Data/results_data_Flevoland_100_20.Rdata")

