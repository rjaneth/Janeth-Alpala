# 
 rm(list = ls())
# if(!require("rstudioapi")) install("rstudioapi")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# for normal distribution
library(ggplot2)
library(e1071)
library(nortest)

#load("./Data/results_data_Ottawa_100_5.Rdata") # difference_value
#load("./Data/results_data_Flevoland_300_5.Rdata")

load("./Data/results_data_Flevoland_300_7.Rdata")
#load("./Data/results_data_Flevoland_300_9.Rdata")

calculate_p_values_matrix <- function(data_matrix, sigma) {
  rows <- nrow(data_matrix)
  cols <- ncol(data_matrix)
  
  p_values_matrix <- matrix(NA, nrow = rows, ncol = cols)
  
  for (i in 1:rows) {
    for (j in 1:cols) {
      test_difference <- data_matrix[i, j]
      
      epsilon <- (test_difference - 0.0) / sigma
      p_value <-2 * pnorm(-abs(epsilon))  
      
      p_values_matrix[i, j] <- p_value
    }
  }
  
  return(p_values_matrix )
}


mean_difference_values <- mean(difference_values, na.rm = TRUE)
sd_difference_values <- sd(difference_values, na.rm = TRUE)


p_values_matrix <- calculate_p_values_matrix(difference_values, sd_difference_values)

save(p_values_matrix, file = "./Data/results_pvalue_Flevoland_300_7_n.Rdata")

aspect_ratio <- ncol(p_values_matrix) / nrow(p_values_matrix)

# Create a sequence of values for x and y axes based on the aspect ratio
x_values <- seq(0, 1, length.out = ncol(p_values_matrix) + 1)
y_values <- seq(0, aspect_ratio, length.out = nrow(p_values_matrix) + 1)

# Plot the p_values_matrix with adjusted axes
image(x = x_values, y = y_values, z = p_values_matrix, 
      col = colorRampPalette(c("navy", "yellow", "orange"))(100),
      main = "Heatmap de p-values",
      xlab = "Columnas",
      ylab = "Filas")

# Add grid lines to visualize the pixel positions
abline(v = x_values, h = y_values, col = "gray", lty = 2)

# Calculate the aspect ratio of the image
# aspect_ratio <- ncol(p_values_matrix) / nrow(p_values_matrix)
# 
# # Invert the y-axis values
# y_values <- rev(seq(0, aspect_ratio, length.out = nrow(p_values_matrix) + 1))
# 
# # Plot the p_values_matrix with adjusted axes
# image(x = seq(0, 1, length.out = ncol(p_values_matrix) + 1), y = y_values, z = p_values_matrix,
#       col = colorRampPalette(c("navy", "yellow", "red"))(100),
#       main = "Heatmap de p-values",
#       xlab = "Columnas",
#       ylab = "Filas")
# aspect_ratio <- ncol(p_values_matrix) / nrow(p_values_matrix)
# 
# # Create a sequence of values for x and y axes based on the aspect ratio
 # x_values <- seq(0, 1, length.out = ncol(p_values_matrix) + 1)
 # # Invertir el eje y al crear y_values
 # # Mantén la definición original de y_values
 # y_values <- seq(0, aspect_ratio, length.out = nrow(p_values_matrix) + 1)
 # 
 # # Invierte el eje y solo para la visualización
 # image(x = x_values, y = rev(y_values), z = p_values_matrix,
 #       col = colorRampPalette(c("navy", "yellow", "red"))(100),
 #       main = "Heatmap de p-values",
 #       xlab = "Columnas",
 #       ylab = "Filas")

# Mantén la definición original de x_values e y_values
# Mantén la definición original de x_values e y_values
# Mantén la definición original de x_values e y_values

# library(gplots)
# library(fields)
# 
# 
# # Mantén la definición original de x_values e y_values
# x_values <- seq(0, 1, length.out = ncol(p_values_matrix) + 1)
# y_values <- seq(0, aspect_ratio, length.out = nrow(p_values_matrix) + 1)
# 
# # Establece parámetros de márgenes para dejar espacio para la leyenda
# par(mar = c(5, 4, 4, 7))
# 
# # Invierte el eje y y el eje z solo para la visualización
# image(x = x_values, y = y_values, z = p_values_matrix[, nrow(p_values_matrix):1],
#       col = colorRampPalette(c("navy", "yellow", "red"))(100),
#       main = "Heatmap de p-values",
#       xlab = "Columnas",
#       ylab = "Filas",
#       axes = FALSE)  # Desactiva los ejes originales
# 
# # Agrega los ejes invertidos
# axis(2, at = seq(aspect_ratio, 0, length.out = nrow(p_values_matrix) + 1), 
#      labels = rev(seq(aspect_ratio, 0, length.out = nrow(p_values_matrix) + 1)), las = 1)
# axis(1)
# 
# # Restaura los parámetros de márgenes a su valor original
# par(mar = c(5, 4, 4, 2) + 0.1)
# 
# # Agrega la barra de colores en gradiente como leyenda, ajustando posición y tamaño
# image.plot(zlim = c(0, 1), legend.only = TRUE, 
#            col = colorRampPalette(c("navy", "yellow", "red"))(100),
#            legend.mar = 0.1, legend.width = 0.1)
# 












#  y_values <- seq(0, aspect_ratio, length.out = nrow(p_values_matrix) + 1)
# # 
# # Plot the p_values_matrix with adjusted axes
# image(x = x_values, y = y_values, z = p_values_matrix,
#       col = colorRampPalette(c("navy", "yellow", "red"))(100),
#       main = "Heatmap de p-values",
#       xlab = "Columnas",
#       ylab = "Filas")

# Add grid lines to visualize the pixel positions
#abline(v = x_values, h = y_values, col = "white", lty = 2)
# image(p_values_matrix, 
#       col = colorRampPalette(c("navy", "yellow", "orange"))(100),
#       main = "Heatmap de p-values",
#       xlab = "Columnas",
#       ylab = "Filas")



# image(p_values_matrix, 
#       col = colorRampPalette(c("blue", "orange"))(100),
#       main = "Heatmap de p-values",
#       xlab = "Columnas",
#       ylab = "Filas")

#source("../imagematrix.R")



#hist(p_values_matrix)
#par(mfrow=c(1,2))
#plot(imagematrix(equalize(x)))
#plot(imagematrix(p_values_matrix>.01))

#imagematrixPNG(imagematrix(equalize(x)), name = "Intensity.png")
#imagematrixPNG(imagematrix(p_values_matrix>0.01), name="pvalue7x7-5.png")








# # Convertir la matriz de p-values a un dataframe para ggplot2
# p_values_df <- as.data.frame(as.table(p_values_matrix))
# colnames(p_values_df) <- c("Row", "Column", "P_Value")
# 
# # Crear el grC!fico con ggplot2
# ggplot(p_values_df, aes(x = as.factor(Column), y = as.factor(Row), fill = P_Value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "blue") +
#   labs(title = "Matrix de P-Values", x = "Columna", y = "Fila") +
#   scale_x_discrete(name = "Columna", labels = seq(1, ncol(p_values_matrix))) +
#   scale_y_discrete(name = "Fila", labels = seq(1, nrow(p_values_matrix))) +
#   theme_minimal()


# p_values_df <- as.data.frame(as.table(p_values_matrix))
# colnames(p_values_df) <- c("Row", "Column", "P_Value")
# 
# # Crear el grC!fico con ggplot2
# ggplot(p_values_df, aes(x = Column, y = Row, fill = P_Value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "blue") +
#   labs(title = "Matrix de P-Values", x = "Columna", y = "Fila") +
#   theme_minimal()