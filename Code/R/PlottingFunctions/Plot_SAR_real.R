
rm(list = ls())

library(raster)
if(!require("rstudioapi")) install("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../../Code/R/Programs/read_ENVI_images.R")


x <- myread.ENVI(file='../../../Data/SAR/Ottawa/Intensity_VV.img', headerfile='../../../Data/SAR/Ottawa/Intensity_VV.hdr')

nCols_Ottawa <- 512
nRows_Ottawa <- 512

centerRow_Ottawa <- floor(nRows_Ottawa / 2)
centerCol_Ottawa <- floor(nCols_Ottawa / 2)

roiSize_Ottawa <- 200

roiUpperStartRow <- 1
roiUpperEndRow <- roiUpperStartRow + roiSize_Ottawa - 1

roiUpperStartCol <- centerCol_Ottawa - floor(roiSize_Ottawa / 2) + 1
roiUpperEndCol <- roiUpperStartCol + roiSize_Ottawa - 1

# ROI en la parte superior
roi <- x[roiUpperStartRow:roiUpperEndRow, roiUpperStartCol:roiUpperEndCol]
escala_log <- 0.01
subx_log <- log10(roi + escala_log)

# Crear un objeto raster
raster_obj <- raster(subx_log)
media_valor2 <- mean(raster_obj)

# Guardar el ROI como PNG
savepng <- function(imagen, nombre) {
  png(file = nombre, width = dim(imagen)[2], height = dim(imagen)[1])
  vectorzeros <- replicate(4, 0)
  par(mar = vectorzeros, oma = vectorzeros, omi = vectorzeros)
  image(imagen, col = gray(seq(0, 1, length.out = 256)))
  #image(imagen, col = gray(seq(0, 0.6, length.out = 51)))
  dev.off()
  return(0)
}

savepng(raster_obj, "./ROIs/Ottawa_roi3.png")

# Crear una copia de la imagen completa
x_copy <- x

# Dibujar un cuadro de contorno azul alrededor del ROI
max_roi_value <- max(roi)

# Dibujar un cuadro de contorno azul alrededor del ROI
x_copy[roiUpperStartRow, roiUpperStartCol:roiUpperEndCol] <- max_roi_value + 1
x_copy[roiUpperEndRow, roiUpperStartCol:roiUpperEndCol] <- max_roi_value + 1
x_copy[roiUpperStartRow:roiUpperEndRow, roiUpperStartCol] <- max_roi_value + 1
x_copy[roiUpperStartRow:roiUpperEndRow, roiUpperEndCol] <- max_roi_value + 1

# Convertir la imagen a escala logarítmica
escala_log_copy <- log10(x_copy + escala_log)

# Crear un objeto raster
raster_obj_copy <- raster(escala_log_copy)

# Ajustar las coordenadas del contorno al tamaño del píxel
pixelSize <- 1  # Tamaño del píxel, ajusta si es necesario

savepng <- function(imagen, nombre) {
  png(file = nombre, width = dim(imagen)[2], height = dim(imagen)[1])
  vectorzeros <- replicate(4, 0)
  par(mar = vectorzeros, oma = vectorzeros, omi = vectorzeros)
  #image(imagen, col = rainbow(51))
  #image(imagen, col = gray(seq(0, 0.6, length.out = 51)))
  image(imagen, col = gray(seq(0, 0.6, length.out = 256)))
  # Dibujar el contorno azul
  rect(roiUpperStartCol - 0.5, roiUpperStartRow - 0.5, roiUpperEndCol + 0.5, roiUpperEndRow + 0.5, border = "blue", lwd = 5)
  dev.off()
  return(0)
}

# Guardar la imagen con el cuadro de contorno azul
savepng(raster_obj_copy, "./ROIs/Ottawa_roi_with_contour3.png")



# x <- myread.ENVI(file='../../../Data/SAR/Ottawa/Intensity_VV.img', headerfile='../../../Data/SAR/Ottawa/Intensity_VV.hdr')
# 
# nCols_Ottawa <- 512
# nRows_Ottawa <- 512
# 
# centerRow_Ottawa <- floor(nRows_Ottawa / 2)
# centerCol_Ottawa <- floor(nCols_Ottawa / 2)
# 
# roiSize_Ottawa <- 50
# 
# roiUpperStartRow <- 1
# roiUpperEndRow <- roiUpperStartRow + roiSize_Ottawa - 1
# 
# roiUpperStartCol <- centerCol_Ottawa - floor(roiSize_Ottawa / 2) + 1
# roiUpperEndCol <- roiUpperStartCol + roiSize_Ottawa - 1
# 
# # ROI en la parte superior
# roi <- x[roiUpperStartRow:roiUpperEndRow, roiUpperStartCol:roiUpperEndCol]
# escala_log <- 0.1
# subx_log <- log10(roi + escala_log)
# 
# # Crear un objeto raster
# raster_obj <- raster(subx_log)
# media_valor2 <- mean(raster_obj)
# 
# # Guardar el ROI como PNG
# savepng <- function(imagen, nombre) {
#   png(file = nombre, width = dim(imagen)[2], height = dim(imagen)[1])
#   vectorzeros <- replicate(4, 0)
#   par(mar = vectorzeros, oma = vectorzeros, omi = vectorzeros)
#   image(imagen, col = gray(seq(0, 0.6, length.out = 51)))
#   dev.off()
#   return(0)
# }
# 
# savepng(raster_obj, "./ROIs/Ottawa_roi1.png")
# 
# # Crear una copia de la imagen completa
# x_copy <- x
# 
# # Dibujar un cuadro de contorno rojo alrededor del ROI
# x_copy[roiUpperStartRow, roiUpperStartCol:roiUpperEndCol] <- max(x_copy) + 1  # Valor fuera del rango de intensidad de la imagen
# x_copy[roiUpperEndRow, roiUpperStartCol:roiUpperEndCol] <- max(x_copy) + 1
# x_copy[roiUpperStartRow:roiUpperEndRow, roiUpperStartCol] <- max(x_copy) + 1
# x_copy[roiUpperStartRow:roiUpperEndRow, roiUpperEndCol] <- max(x_copy) + 1
# 
# # Convertir la imagen a escala logarítmica
# escala_log_copy <- log10(x_copy + escala_log)
# 
# # Crear un objeto raster
# raster_obj_copy <- raster(escala_log_copy)
# 
# # Ajustar las coordenadas del contorno al tamaño del píxel
# pixelSize <- 1  # Tamaño del píxel, ajusta si es necesario
# 
# # Guardar la imagen con el cuadro de contorno rojo
# savepng(raster_obj_copy, "./ROIs/Ottawa_roi_with_contour1.png")



# library(raster)
# 
# source("../../../Code/R/Programs/read_ENVI_images.R")
# 
# 
# 
# x <- myread.ENVI(file='../../../Data/SAR/Ottawa/Intensity_VV.img', headerfile='../../../Data/SAR/Ottawa/Intensity_VV.hdr')
# 
# nCols_Ottawa <- 512
# nRows_Ottawa <- 512
# 
# centerRow_Ottawa <- floor(nRows_Ottawa / 2)
# centerCol_Ottawa <- floor(nCols_Ottawa / 2)
# 
# roiSize_Ottawa <- 200
# 
# roiUpperStartRow <- 1
# roiUpperEndRow <- roiUpperStartRow + roiSize_Ottawa - 1
# 
# roiUpperStartCol <- centerCol_Ottawa - floor(roiSize_Ottawa / 2) + 1
# roiUpperEndCol <- roiUpperStartCol + roiSize_Ottawa - 1
# 
# # ROI en la parte superior
# roi <- x[roiUpperStartRow:roiUpperEndRow, roiUpperStartCol:roiUpperEndCol]
# escala_log <- 0.1
# subx_log <- log10(roi + escala_log)
# 
# # Crear un objeto raster
# raster_obj <- raster(subx_log)
# media_valor2 <- mean(raster_obj)
# 
# # Guardar el ROI como PNG
# savepng <- function(imagen, nombre) {
#   png(file = nombre, width = dim(imagen)[2], height = dim(imagen)[1])
#   vectorzeros <- replicate(4, 0)
#   par(mar = vectorzeros, oma = vectorzeros, omi = vectorzeros)
#   image(imagen, col = gray(seq(0, 0.6, length.out = 51)))
#   dev.off()
#   return(0)
# }
# 
# savepng(raster_obj, "./ROIs/Ottawa_roi.png")
# 
# # Crear una copia de la imagen completa
# x_copy <- x
# 
# # Dibujar un cuadro de contorno rojo alrededor del ROI
# x_copy[roiUpperStartRow, roiUpperStartCol:roiUpperEndCol] <- 1
# x_copy[roiUpperEndRow, roiUpperStartCol:roiUpperEndCol] <- 1
# x_copy[roiUpperStartRow:roiUpperEndRow, roiUpperStartCol] <- 1
# x_copy[roiUpperStartRow:roiUpperEndRow, roiUpperEndCol] <- 1
# 
# # Convertir la imagen a escala logarítmica
# escala_log_copy <- log10(x_copy + escala_log)
# 
# # Crear un objeto raster
# raster_obj_copy <- raster(escala_log_copy)
# 
# # Guardar la imagen con el cuadro de contorno rojo
# # savepng <- function(imagen, nombre, roiStartRow, roiStartCol, roiEndRow, roiEndCol) {
# #   png(file = nombre, width = dim(imagen)[2], height = dim(imagen)[1])
# #   vectorzeros <- replicate(4, 0)
# #   par(mar = vectorzeros, oma = vectorzeros, omi = vectorzeros)
# #   
# #   # Cambiar la paleta de colores a un espectro de colores, como "rainbow"
# #   image(imagen, col = rainbow(51))
# #   
# #   
# #   # Dibujar el contorno rojo
# #   rect(roiStartCol - 0.5, roiStartRow - 0.5, roiEndCol + 0.5, roiEndRow + 0.5, border = "red", lwd = 2)
# #   
# #   dev.off()
# #   return(0)
# # }
# savepng <- function(imagen, nombre, roiStartRow, roiStartCol, roiEndRow, roiEndCol) {
#   png(file = nombre, width = dim(imagen)[2], height = dim(imagen)[1])
#   vectorzeros <- replicate(4, 0)
#   par(mar = vectorzeros, oma = vectorzeros, omi = vectorzeros)
#   
#   # Dibujar la imagen en escala de grises
#   image(imagen, col = gray(seq(0, 0.6, length.out = 51)))
#   
#   # Ajustar las coordenadas del contorno al tamaño del píxel
#   pixelSize <- 1  # Tamaño del píxel, ajusta si es necesario
#   
#   # Dibujar el contorno rojo
#   rect(roiStartCol - 0.5 * pixelSize, roiStartRow - 0.5 * pixelSize,
#        roiEndCol + 0.5 * pixelSize, roiEndRow + 0.5 * pixelSize,
#        border = "blue", lwd = 2)
#   
#   dev.off()
#   return(0)
# }
# 
# 
# 
# 
# # Guardar la imagen con el cuadro de contorno rojo
# savepng(raster_obj_copy, "./ROIs/Ottawa_roi_with_contour.png", roiUpperStartRow, roiUpperStartCol, roiUpperEndRow, roiUpperEndCol)
# 






# x <- myread.ENVI(file='../../../Data/SAR/Ottawa/Intensity_VV.img', headerfile='../../../Data/SAR/Ottawa/Intensity_VV.hdr')
# 
# nCols_Ottawa <- 512
# nRows_Ottawa <- 512
# 
# centerRow_Ottawa <- floor(nRows_Ottawa / 2)
# centerCol_Ottawa <- floor(nCols_Ottawa / 2)
# 
# roiSize_Ottawa <- 100
# 
# roiUpperStartRow <- 1
# roiUpperEndRow <- roiUpperStartRow + roiSize_Ottawa - 1
# 
# roiUpperStartCol <- centerCol_Ottawa - floor(roiSize_Ottawa / 2) + 1
# roiUpperEndCol <- roiUpperStartCol + roiSize_Ottawa - 1
# 
# # ROI en la parte superior
# roi <- x[roiUpperStartRow:roiUpperEndRow, roiUpperStartCol:roiUpperEndCol]
# escala_log <- 0.1
# subx_log <- log10(roi + escala_log)
# 
# # Crear un objeto raster
# raster_obj <- raster(subx_log)
# media_valor2 <- mean(raster_obj)
# 
# # Guardar el ROI como PNG
# savepng <- function(imagen, nombre) {
#   png(file = nombre, width = dim(imagen)[2], height = dim(imagen)[1])
#   vectorzeros <- replicate(4, 0)
#   par(mar = vectorzeros, oma = vectorzeros, omi = vectorzeros)
#   image(imagen, col = gray(seq(0, 0.6, length.out = 51)))
#   dev.off()
#   return(0)
# }
# 
# savepng(raster_obj, "./ROIs/Ottawa_roi.png")
# 
# # Crear una copia de la imagen completa
# x_copy <- x
# 
# # Dibujar un cuadro de contorno rojo alrededor del ROI
# x_copy[roiUpperStartRow, roiUpperStartCol:roiUpperEndCol] <- 1
# x_copy[roiUpperEndRow, roiUpperStartCol:roiUpperEndCol] <- 1
# x_copy[roiUpperStartRow:roiUpperEndRow, roiUpperStartCol] <- 1
# x_copy[roiUpperStartRow:roiUpperEndRow, roiUpperEndCol] <- 1
# 
# # Convertir la imagen a escala logarítmica
# escala_log_copy <- log10(x_copy + escala_log)
# 
# # Crear un objeto raster
# raster_obj_copy <- raster(escala_log_copy)
# 
# # Guardar la imagen con el cuadro de contorno rojo
# savepng <- function(imagen, nombre) {
#   png(file = nombre, width = dim(imagen)[2], height = dim(imagen)[1])
#   vectorzeros <- replicate(4, 0)
#   par(mar = vectorzeros, oma = vectorzeros, omi = vectorzeros)
#   image(imagen, col = gray(seq(0, 0.6, length.out = 51)))
#   # Dibujar el contorno rojo
#   rect(roiUpperStartCol - 0.5, roiUpperStartRow - 0.5, roiUpperEndCol + 0.5, roiUpperEndRow + 0.5, border = "red", lwd = 2)
#   dev.off()
#   return(0)
# }
# 
# # Guardar la imagen con el cuadro de contorno rojo
# savepng(raster_obj_copy, "./ROIs/Ottawa_roi_with_contour.png")




# #x_Brazil <- myread.ENVI(file='../../../Data/SAR/Brazil/Intensity_VV.img', headerfile='../../../Data/SAR/Brazil/Intensity_VV.hdr')
# #x_Berlin<- myread.ENVI(file='../../../Data/SAR/Berlin/Intensity_VV.img', headerfile='../../../Data/SAR/Berlin/Intensity_VV.hdr')
# #x_Flevoland <- myread.ENVI(file='../../../Data/SAR/Flevoland/Intensity_VV.img', headerfile='../../../Data/SAR/Flevoland/Intensity_VV.hdr')
# x <- myread.ENVI(file='../../../Data/SAR/Ottawa/Intensity_VV.img', headerfile='../../../Data/SAR/Ottawa/Intensity_VV.hdr')
# #x <- myread.ENVI(file='./Ottawa/Intensity_VV.img', headerfile='./Ottawa/Intensity_VV.hdr')
# 
# 
# nCols_Ottawa = 512
# nRows_Ottawa= 512
# 
# centerRow_Ottawa  <- floor(nRows_Ottawa / 2)
# centerCol_Ottawa  <- floor(nCols_Ottawa  / 2)
# 
# roiSize_Ottawa <- 512
# # 
# 
# 
# roiUpperStartRow <- 1
# roiUpperEndRow <- roiUpperStartRow + roiSize_Ottawa - 1
# # 
# roiUpperStartCol <- centerCol_Ottawa - floor(roiSize_Ottawa / 2) + 1
# roiUpperEndCol <- roiUpperStartCol + roiSize_Ottawa - 1
# # 
# #  ROI en la parte superior
# roi<- x[roiUpperStartRow:roiUpperEndRow, roiUpperStartCol:roiUpperEndCol]
# escala_log <- 0.1
# subx_log <- log10(roi + escala_log)
# 
# # Crear un objeto raster
# 
# raster_obj <- raster(subx_log)
# media_valor2 <- mean(raster_obj)
# 
# 
# savepng <- function(imagen, nombre) {
#   png(file = nombre, width = dim(imagen)[2], height = dim(imagen)[1])
#   vectorzeros <- replicate(4, 0)
#   par(mar = vectorzeros, oma = vectorzeros, omi = vectorzeros)
#   image(imagen, col = gray(seq(0, 0.6, length.out = 51)))
#   dev.off()
#   return(0)
# }
# 
# # 
# savepng(raster_obj, "./ROIs/Ottawa_roi.png")
# 
