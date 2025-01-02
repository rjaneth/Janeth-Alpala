# 408
# Inicia el temporizador
start_time <- Sys.time()

# Tu código aquí
source("../../../Code/R/MainFunctions/bootstrap_al_omari_1_estimator.R")
source("../../../Code/R/MainFunctions/bootstrap_correa_estimator_log_mean.R")
source("../../../Code/R/Programs/read_ENVI_images.R")
source("../../../Code/R/MainFunctions/al_omari_1_estimator.R")
source("../../../Code/R/MainFunctions/correa_estimator.R")

# x <- myread.ENVI(file='../../../Data/SAR/Ottawa_512/Intensity_VV.img',
#                  headerfile='../../../Data/SAR/Ottawa_512/Intensity_VV.hdr')
 # x <- myread.ENVI(file='../../../Data/SAR/mexico_512/Intensity_VV.img', 
 #                  headerfile='../../../Data/SAR/mexico_512/Intensity_VV.hdr')
 
 # x <- myread.ENVI(file='../../../Data/SAR/chicago_1024/Intensity_HH.img', 
 #                  headerfile='../../../Data/SAR/chicago_1024/Intensity_HH.hdr')
 # 
 # x <- myread.ENVI(file='../../../Data/SAR/Ottawa_1024/Intensity_VV.img', 
 #                  headerfile='../../../Data/SAR/Ottawa_1024/Intensity_VV.hdr')
# x <- myread.ENVI(file='../../../Data/SAR/Frankfurt_1024/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/Frankfurt_1024/Intensity_VV.hdr')
# x <- myread.ENVI(file='../../../Data/SAR/mexico600/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/mexico600/Intensity_VV.hdr')
# x <- myread.ENVI(file='../../../Data/SAR/chicago_urban_1024/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/chicago_urban_1024/Intensity_VV.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/SanFrancisco_650/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/SanFrancisco_650/Intensity_VV.hdr')
# x <- myread.ENVI(file='../../../Data/SAR/Chicago_400/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/Chicago_400/Intensity_VV.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/Lake_512/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/Lake_512/Intensity_VV.hdr')

 x <- myread.ENVI(file='../../../Data/SAR/panama512/Intensity_VV.img', 
                   headerfile='../../../Data/SAR/panama512/Intensity_VV.hdr')

#load("../Programs/Data/Phantom_4_regions.Rdata")

# Asignar las dimensiones de la matriz cargada a las variables rows y cols
#save(Z, file="./Data/Phantom_4_regions.Rdata")

# rows <- nrow(Z)
# cols <- ncol(Z)


#C:/Users/luiso/Documents/Github/Janeth-Alpala/Data/SAR/chicago_urban_1024

L <- 5
B <- 200
window_size <- 7

rows <- nrow(x)
cols <- ncol(x)
difference_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)

# Iterar sobre ventanas deslizantes
for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Seleccionar ventana local
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
    difference_values[i, j] <- bootstrap_al_omari_1_estimator(window_data,B) - (log(mean(window_data)) + (L - log(L) + lgamma(L) + (1 - L) * digamma(L)))
    #difference_values[i, j] <- correa_estimator(window_data) 
  }
}

# Guarda los resultados
save(difference_values, x, file = "./Data/results_Panama_512_AO_7x7_L5_200b.Rdata")

# Detén el temporizador
end_time <- Sys.time()

# Calcula el tiempo transcurrido
execution_time <- end_time - start_time
execution_time

