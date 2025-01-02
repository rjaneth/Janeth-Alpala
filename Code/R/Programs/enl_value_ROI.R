
rm(list = ls())
# Cargar la imagen SAR

source("../../../Code/R/Programs/read_ENVI_images.R")

# Cargar la imagen SAR
x <- myread.ENVI(file='../../../Data/SAR/Houston_1024/Intensity_HH.img', 
                 headerfile='../../../Data/SAR/Houston_1024/Intensity_HH.hdr')

# Definir las coordenadas de la región de interés (ROI)
# Ejemplo: una región homogénea en el centro de la imagen
roi_start_row <- floor(nrow(x) / 2) - 3  # Centrado, con ventana de 7x7
roi_start_col <- floor(ncol(x) / 2) - 3  # Centrado, con ventana de 7x7

# Seleccionar la ventana local
window_data <- x[roi_start_row:(roi_start_row + 6), roi_start_col:(roi_start_col + 6)]

# Calcular la media muestral y la varianza muestral en la ventana
mean_value <- mean(window_data)
variance_value <- var(window_data)

# Calcular el ENL en la ventana
enl_value <- (mean_value^2) / variance_value

# Imprimir el valor de ENL
print(paste("ENL for the selected ROI: ", enl_value))

# Guardar los valores de ENL y la imagen original
save(enl_value, x, file = "./Data/ENL_results_Houston_1024_ROI.Rdata")

# Cargar los valores de ENL y la imagen original desde el archivo guardado
load("./Data/ENL_results_Houston_1024_ROI.Rdata")

# Imprimir el valor de ENL cargado
print(enl_value)