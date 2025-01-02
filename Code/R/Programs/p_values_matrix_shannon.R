#556 
#rm(list = ls())

load("./Data/results_Phantom_4_AO_B100_w7_L5.Rdata")
calculate_p_values_matrix <- function(data_matrix, mu, sigma) {
  rows <- nrow(data_matrix)
  cols <- ncol(data_matrix)
  
  p_values_matrix <- matrix(NA, nrow = rows, ncol = cols)
  
  for (i in 1:rows) {
    for (j in 1:cols) {
      test_difference <- data_matrix[i, j]
      
      epsilon <-(test_difference - 0) / (sigma)
      p_value <-  2*pnorm(-abs(epsilon))
  
      p_values_matrix[i, j] <- p_value
    }
  }
  
  return(p_values_matrix )
}




mean_difference_values <- mean(difference_values_shannon, na.rm = TRUE)
sd_difference_values <- sd(difference_values_shannon, na.rm = TRUE)


p_values_matrix <- calculate_p_values_matrix(difference_values_shannon, mean_difference_values, sd_difference_values)

#save(p_values_matrix, file = "./Data/results_pvalue_Phantom_4_regions_AO_7W_L5_100b.Rdata")

source("../imagematrix.R")

#plot(imagematrix(p_values_matrix ))
plot(imagematrix(p_values_matrix> 0.05))
# plot(imagematrix(equalize(difference_values_shannon)))


# 
# #     imagematrixPNG(imagematrix(equalize(x)), name = "sinaloa_512.png")
imagematrixPNG(imagematrix(equalize(difference_values)), name = "Entropy_Illinois_2_renyi_1024_w7_b100_L36_098.png")
# # # # # # # # # # # 
imagematrixPNG(imagematrix(p_values_matrix), name="H_pvalue_Illinois_2_renyi_1024_w7_b100_L36_098.png")
imagematrixPNG(imagematrix(p_values_matrix>0.05), name="H_005_Illinois_2_renyi_1024_w7_b100_L36_098.png")
