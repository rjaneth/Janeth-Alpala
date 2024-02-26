rm(list = ls())
# if(!require("rstudioapi")) install("rstudioapi")
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# for normal distribution
library(ggplot2)
library(e1071)
library(nortest)
library(fields)
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
      
      epsilon <-(test_difference - 0.0) / (sigma)
      p_value <-p_value <- 2 * pnorm(abs(epsilon))-1#  
      
      p_values_matrix[i, j] <- p_value
    }
  }
  
  return(p_values_matrix )
}




mean_difference_values <- mean(difference_values, na.rm = TRUE)
sd_difference_values <- sd(difference_values, na.rm = TRUE)


p_values_matrix <- calculate_p_values_matrix(difference_values, sd_difference_values)


source("../imagematrix.R")
#hist(p_values_matrix)
plot(imagematrix(p_values_matrix),significance_level = 0.05)

