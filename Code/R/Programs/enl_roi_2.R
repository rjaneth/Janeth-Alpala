rm(list = ls())

source("../../../Code/R/Programs/read_ENVI_images.R")

# Cargar la imagen SAR
x <- myread.ENVI(file='../../../Data/SAR/Houston_1024/Intensity_HH.img', 
                 headerfile='../../../Data/SAR/Houston_1024/Intensity_HH.hdr')

# Definir el tamaño de la ventana
window_size <- 7

# Obtener las dimensiones de la imagen
nRows <- nrow(x)
nCols <- ncol(x)

# Definir las coordenadas del centro de la imagen
centerRow <- floor(nRows / 2)
centerCol <- floor(nCols / 2)

# Calcular las coordenadas del ROI (ventana central)
roiStartRow <- centerRow - floor(window_size / 2) + 1
roiEndRow <- roiStartRow + window_size - 1
roiStartCol <- centerCol - floor(window_size / 2) + 1
roiEndCol <- roiStartCol + window_size - 1

# Extraer el ROI
roi <- x[roiStartRow:roiEndRow, roiStartCol:roiEndCol]

# Calcular la media y la varianza en la ventana
mean_value <- mean(roi)
variance_value <- var(as.vector(roi))

# Calcular el ENL en la ventana
enl_value <- (mean_value^2) / variance_value

# Verificar si enl_value es un solo valor y corregir si es necesario
# if (length(enl_value) != 1) {
#   stop("Error in ENL calculation: ENL value is not a single number.")
# }

# Imprimir el valor de ENL
print(paste("ENL for the selected ROI: ", enl_value))

# Guardar los valores de ENL y la imagen original
save(enl_value, x, file = "./Data/ENL_results_Houston_1024_Central_ROI.Rdata")

# Cargar los valores de ENL y la imagen original desde el archivo guardado
load("./Data/ENL_results_Houston_1024_Central_ROI.Rdata")

# Imprimir el valor de ENL cargado
print(enl_value)
