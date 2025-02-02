source("../../../Code/R/Programs/read_ENVI_images.R")

# Cargar la imagen SAR
# x <- myread.ENVI(file='../../../Data/SAR/agua_300/Intensity_VV.img',
#                  headerfile='../../../Data/SAR/agua_300/Intensity_VV.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/Houston_100/Intensity_HH.img', 
#                  headerfile='../../../Data/SAR/Houston_100/Intensity_HH.hdr')

# x <- myread.ENVI(file='../../../Data/SAR/envi_sf_water_200/Intensity_VV.img', 
#                  headerfile='../../../Data/SAR/envi_sf_water_200/Intensity_VV.hdr')

#x <- myread.ENVI(file='../../../Data/SAR/MUNICH_80/Intensity_HH.img',
#                 headerfile='../../../Data/SAR/MUNICH_80/Intensity_HH.hdr')
x <- myread.ENVI(file='../../../Data/SAR/sinaloa_512/Intensity_VV.img', 
                 headerfile='../../../Data/SAR/sinaloa_512/Intensity_VV.hdr')

window_size <-20


rows <- nrow(x)
cols <- ncol(x)


enl_values <- matrix(NA, nrow = rows - window_size + 1, ncol = cols - window_size + 1)


for (i in 1:(rows - window_size + 1)) {
  for (j in 1:(cols - window_size + 1)) {
    # Seleccionar la ventana local
    window_data <- x[i:(i + window_size - 1), j:(j + window_size - 1)]
    
    
    mean_value <- mean(as.vector(window_data))
    variance_value <- var(as.vector(window_data))
    
    
    
    
    enl_value <- (mean_value^2) / variance_value
    
    # Verificar si enl_value es un solo valor
    if (length(enl_value) == 1) {
      # Almacenar el valor de ENL en la matriz
      enl_values[i, j] <- enl_value
    }
  }
}


#save(enl_values, x, file = "./Data/ENL_results_SF_11.Rdata")

# Cargar los valores de ENL y la imagen original desde el archivo guardado
#load("./Data/ENL_results_SF_11.Rdata")

# Calcular el promedio de todos los valores de ENL
mean_enl <- mean(enl_value, na.rm = TRUE)

# Imprimir el promedio de ENL
print(mean_enl)

