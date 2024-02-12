
rm(list = ls())
if(!require("rstudioapi")) install("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
      
      epsilon <- test_difference / sigma
      p_value <- 2 * pnorm(abs(epsilon)) - 1
      
      p_values_matrix[i, j] <- p_value
    }
  }
  
  return(p_values_matrix)
}


mean_difference_values <- mean(difference_values, na.rm = TRUE)
sd_difference_values <- sd(difference_values, na.rm = TRUE)


p_values_matrix <- calculate_p_values_matrix(difference_values, sd_difference_values)


# p_values_df <- as.data.frame(as.table(p_values_matrix))
# colnames(p_values_df) <- c("Row", "Column", "P_Value")
# 
# 
# ggplot(p_values_df, aes(x = as.factor(Column), y = as.factor(Row), fill = P_Value)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "blue") +
#   labs(title = "Matrix de P-Values") +
#   theme_minimal() +
#   theme(axis.text = element_blank(), axis.title = element_blank())

source("../imagematrix.R")

hist(p_values_matrix)
par(mfrow=c(1,2))
plot(imagematrix(equalize(x)))
plot(imagematrix(p_values_matrix))

imagematrixPNG(imagematrix(equalize(x)), name = "Intensity.png")
imagematrixPNG(imagematrix(p_values_matrix), name="pvalue7x7.png")




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