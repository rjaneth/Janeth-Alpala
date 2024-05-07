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
load("./Data/results_Phantom_4_z1_200b.Rdata")#OK
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

#save(p_values_matrix, file = "./Data/results_pvalue_Phantom_4_regions_AO_7W_L5_100b.Rdata")

source("../imagematrix.R")
hist(p_values_matrix)
#plot(imagematrix(p_values_matrix),significance_level = 0.05)
plot(imagematrix(p_values_matrix ))
plot(imagematrix(p_values_matrix >0.1))
plot(imagematrix(equalize(difference_values)))
#plot(imagematrix(equalize(x)))


#   imagematrixPNG(imagematrix(equalize(x)), name = "lake_512.png")
    imagematrixPNG(imagematrix(equalize(difference_values)), name = "Entropy_Phantom_4_z1_200.png")
# # # # 
    imagematrixPNG(imagematrix(p_values_matrix), name="H_pvalue_Phantom_4_z1_200b.png")
    imagematrixPNG(imagematrix(p_values_matrix>0.05), name="H_005_Phantom_4_z1_AO_200b.png")
