rm(list = ls())

#load("./Data/results_Mexico_512_7_AO_300b.Rdata")#
#load("./Data/results_Phantom_4_regions_AO_7W_L5_100b.Rdata")
#load("./Data/results_Chicago_1024_9_AO_200b_36L.Rdata")#OK
#load("./Data/results_Mexico_512_7_AO_300b.Rdata")#OK
#load("./Data/results_panama_512_AO_9W_L5_100b.Rdata")
#save(difference_values, x, file = "./Data/results_Panama_512_AO_7x7_L5_200b.Rdata")
#load("./Data/results_Panama_512_AO_7x7_L5_200b.Rdata")
#load("./Data/results_Lake_512_9W_AO_100b_36L.Rdata")#OK
#save(difference_values, Z, file = "./Data/results_Phantom_4_z1_200b.Rdata")
#load("./Data/results_Phantom_4_z1_200b.Rdata")#OK

#save(difference_values, x, file = "./Data/results_Michigan_1024_7_AO_L1.Rdata")
#load("./Data/results_lake_envi_1_9.Rdata")#OK
#load("./Data/results_Rotterdam_AO_L1_7.Rdata")#OK
#load("./Data/results_Rotterdam_AO_L1_7_2.Rdata")#OK

#save(difference_values, x, file = "./Data/results_Rotterdam_AO_L1_7_withou_B.Rdata")
#load("./Data/results_Rotterdam_AO_L1_7_withou_B.Rdata")#OK

#save(difference_values, x, file = "./Data/results_Rotterdam_ebrahimi_L1_7_withou_B.Rdata")
#load("./Data/results_Rotterdam_ebrahimi_L1_7_withou_B.Rdata")#OK
#save(difference_values, x, file = "./Data/results_Rotterdam_AO2_L1_7n_withou_B.Rdata")
#save(difference_values, mean_difference, sd_difference, file = "./Data/results_agua_envi_11.Rdata")
#load("./Data/results_agua_envi_11.Rdata")#OK
#results_Baja_California_1024_7_AO.Rdata
#save(difference_values, x, file = "./Data/results_mato_grosso_2048_7_AO.Rdata")
load("./Data/results_mato_grosso_1258_7_AO.Rdata")#OK
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




mean_difference_values <- mean(filtered_difference_values, na.rm = TRUE)
sd_difference_values <- sd(filtered_difference_values, na.rm = TRUE)


p_values_matrix <- calculate_p_values_matrix(filtered_difference_values, mean_difference_values, sd_difference_values)

#save(p_values_matrix, file = "./Data/results_pvalue_Phantom_4_regions_AO_7W_L5_100b.Rdata")
# Eliminar filas con NA en p_values_matrix
p_values_matrix <- na.omit(p_values_matrix)
source("../imagematrix.R")
#hist(p_values_matrix)
#plot(imagematrix(p_values_matrix),significance_level = 0.05)
#plot(imagematrix(p_values_matrix ))
plot(imagematrix(p_values_matrix >0.1))
#plot(imagematrix(equalize(difference_values)))
#plot(imagematrix(equalize(x)))


 #imagematrixPNG(imagematrix(equalize(x)), name = "mato_grosso_1258.png")
#         imagematrixPNG(imagematrix(equalize(difference_values)), name = "Entropy_AO_Rotterdam_1024_B_ebrahimi.png")
# # # # # # # # #
#       imagematrixPNG(imagematrix(p_values_matrix), name="H_pvalue_ebrahimi_Rotterdam_1024_B.png")
     #  imagematrixPNG(imagematrix(p_values_matrix>0.05), name="H_005_AO2_Rotterdam_1024_B.png")
