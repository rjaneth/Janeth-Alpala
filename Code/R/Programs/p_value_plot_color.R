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

#load("./Data/results_data_Flevoland_300_7.Rdata")
#save(difference_values,z.up.le , file = "./Data/results_data_simulated_z_up_le_entropy_250_7.Rdata")
#load("./Data/results_data_simulated_z_up_le_entropy_250_7.Rdata")
#save(difference_values, x, file = "./Data/results_data_SanFrancisco_650_7.Rdata")
#load("./Data/results_data_SanFrancisco_650_7.Rdata")
#save(difference_values, x, file = "./Data/results_data_Ottawa_512_7_n.Rdata")
#load("./Data/results_data_Ottawa_512_7_n.Rdata")

#save(difference_values, x, file = "./Data/results_data_Flevoland_600_7_n.Rdata")
#load("./Data/results_data_Flevoland_600_7_n.Rdata")

#save(difference_values, x, file = "./Data/results_data_Frankfurt_512_7_n.Rdata")
#save(difference_values, x, file = "./Data/results_data_Flevoland_512_7_n.Rdata")
#save(difference_values, x, file = "./Data/results_data_Flevoland_512_7_n.Rdata")
#save(difference_values, x, file = "./Data/results_chicago_city_512_7_AO_100b.Rdata")
#save(difference_values, x, file = "./Data/results_data_mexico_512_5_Ebr_50b.Rdata")
#save(difference_values, x, file = "./Data/results_Mexico_512_7_AO_300b.Rdata")
#load("./Data/results_Mexico_512_7_AO_300b.Rdata")#OK
#load("./Data/results_Mexico_512_9_AO_100b_6L.Rdata")
#load("./Data/results_Chicago_1024_9_AO_100b_6L.Rdata")
#load("./Data/results_Chicago_1024_9_AO_200b_36L.Rdata")#OK
# load("./Data/results_Mexico_512_9_AO_100b_6L.Rdata")
#load("./Data/results_Frankfurt_1024_9W_AO_100b_5L.Rdata")
# load("./Data/results_Chicago_1024_9_AO_200b_36L.Rdata")
# load("./Data/results_ottawa_city_512_7_AO_300b.Rdata")#OK
# load("./Data/results_Chicago_1024_9_AO_100b_6L.Rdata")
# load("./Data/results_mexico_600_7W_AO_50b_6L.Rdata")
#load("./Data/results_mexico_600_9W_AO_100b_6L.Rdata")
#load("./Data/results_chicago_urban_1024_9W_AO_100b_6L.Rdata")
# load("./Data/results_SanFrancisco_650_9W_AO_100b_5L.Rdata")
# load("./Data/results_Lake_512_9W_AO_100b_36L.Rdata")#OK
#load("./Data/results_Chicago_1024_9_AO_200b_36L.Rdata")
#save(difference_values, x, file = "./Data/results_SanFrancisco_650_9W_AO_100b_5L.Rdata")
#save(difference_values, x, file = "./Data/results_chicago_urban_1024_9W_AO_100b_6L.Rdata")
#results_data_Ottawa_512_7_n
#save(difference_values, x, file = "./Data/results_Ottawa_1024_9W_AO_100b_5L.Rdata")

#save(difference_values, x, file = "./Data/results_Chicago_1024_9_AO_200b_36L.Rdata")
#save(difference_values, x, file = "./Data/results_Chicago_1024_9_AO_100b_6L.Rdata")
#save(difference_values, x, file = "./Data/results_Mexico_512_9_AO_100b_6L.Rdata")
#load("../Programs/Data/results_data_otawa_512_5_Ebr_50b.Rdata")
#load("../Programs/Data/results_ottawa_city_512_7_AO_100b.Rdata")
#load("../Programs/Data/results_data_otawa_512_5_Ebr_50b.Rdata")
#load("../Programs/Data/results_ottawa_city_512_7_AO_300b.Rdata")
#load("../Programs/Data/results_ottawa_city_512_7_Ebrahimi_100b.Rdata")
#load("../Programs/Data/results_ottawa_city_512_7_Ebrahimi_300b.Rdata")
#save(difference_values, x, file = "./Data/results_ottawa_city_512_7_Ebrahimi_300b.Rdata")
#save(difference_values, x, file = "./Data/results_data_otawa_512_5_Ebr_50b.Rdata")
#save(difference_values, x, file = "./Data/results_data_Ottawa_512_7_n.Rdata")
#load("./Data/results_panama_512_AO_9W_L5_100b.Rdata")#results_mexico_600_9W_AO_100b_6L #OK
#load("./Data/results_Mexico_512_7_AO_300b.Rdata")#results_mexico_600_9W_AO_100b_6L
#load("./Data/results_Phantom_4_regions_AO_7W_L5_100b.Rdata")
#C:/Users/luiso/Documents/Github/Janeth-Alpala/Code/R/Programs/Data/results_data_simulated_z_up_le_entropy_250_7_3normal.Rdata
#save(difference_values, x, file = "./Data/results_Phantom_4_regions_AO_7W_L5_100b.Rdata")
# save(difference_values, x, file = "./Data/results_panama_512_AO_9W_L5_100b.Rdata")
#save(difference_values, x, file = "./Data/results_Panama_700_AO_7W_L5_200b.Rdata")
load("./Data/results_Panama_700_AO_7W_L5_200b.Rdata")
calculate_p_values_matrix <- function(data_matrix, mu, sigma) {
  rows <- nrow(data_matrix)
  cols <- ncol(data_matrix)
  
  p_values_matrix <- matrix(NA, nrow = rows, ncol = cols)
  
  for (i in 1:rows) {
    for (j in 1:cols) {
      test_difference <- data_matrix[i, j]
      
      epsilon <-(test_difference - 0) / (sigma)
     p_value <- 2 * pnorm(-abs(epsilon))#  
      
      p_values_matrix[i, j] <- p_value
    }
  }
  
  return(p_values_matrix )
}




mean_difference_values <- mean(difference_values, na.rm = TRUE)
sd_difference_values <- sd(difference_values, na.rm = TRUE)


p_values_matrix <- calculate_p_values_matrix(difference_values, mean_difference_values, sd_difference_values)


source("../imagematrix.R")
#hist(p_values_matrix)
#plot(imagematrix(p_values_matrix),significance_level = 0.05)
plot(imagematrix(p_values_matrix ))
#plot(imagematrix(p_values_matrix <0.1))
plot(imagematrix(p_values_matrix >0.05))
plot(imagematrix(equalize(difference_values)))
plot(imagematrix(equalize(x)))
# imagematrixPNG(imagematrix(equalize(x)), name = "Intensity_Frankfurt_512_n.png")
# #imagematrixPNG(imagematrix(x), name="SanFrancisco_650.png")
# 
 #imagematrixPNG(imagematrix(p_values_matrix), name="Flev512_pvalue7x7_total_H2.png")
 #imagematrixPNG(imagematrix(p_values_matrix>0.1), name="Frankfurt_512_pvalue7x7_H_01.png")
