#409 
# Inicia el temporizador
start_time <- Sys.time()

source("../../../Code/R/MainFunctions/bootstrap_ebrahimi_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator_log_mean.R")
source("../../../Code/R/Programs/read_ENVI_images.R")
source("../../../Code/R/MainFunctions/ebrahimi_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_al_omari_1_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator_log_mean.R")
source("../../../Code/R/Programs/read_ENVI_images.R")
source("../../../Code/R/MainFunctions/al_omari_1_estimator.R")

# x <- myread.ENVI(file='../../../Data/SAR/Ottawa_512/Intensity_VV.img',
#                  headerfile='../../../Data/SAR/Ottawa_512/Intensity_VV.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/Michigan_1024/Intensity_VV.img',
#                  headerfile='../../../Data/SAR/Michigan_1024/Intensity_VV.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/Houston_100/Intensity_HH.img', 
#                  headerfile='../../../Data/SAR/Houston_100/Intensity_HH.hdr')

x <- myread.ENVI(file='../../../Data/SAR/envi_1/Intensity_VV.img', 
                 headerfile='../../../Data/SAR/envi_1/Intensity_VV.hdr')


L <- 36
B <- 100
window_size <- 7
rows <- nrow(x)
cols <- ncol(x)
difference_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterar sobre ventanas deslizantes
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Seleccionar ventana local
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
    difference_values[i, j] <- bootstrap_al_omari_1_estimator(window_data,B)- (log(mean(window_data)) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L)))
  }
}

# Guarda los resultados
save(difference_values, x, file = "./Data/results_lake_envi_1_7.Rdata")
#save(difference_values, x, file = "./Data/results_Houston_100_11_AO_L1.Rdata")

# Detén el temporizador
end_time <- Sys.time()

# Calcula el tiempo transcurrido
execution_time <- end_time - start_time
execution_time

