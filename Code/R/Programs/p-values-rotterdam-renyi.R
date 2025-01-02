rm(list = ls())
#load("../Programs/Data/Phantom_4_z.Rdata")
#save( cd_values_mnad, x, file = "./Data/MnAD_results_data_Rotterdam_1024.Rdata")
#load("../Programs/Data/MnAD_results_data_Rotterdam_1024.Rdata")
#load("../Programs/Data/Phantom_4_regions.Rdata")
#load("./Data/results_Mexico_512_7_AO_300b.Rdata")#
#load("./Data/results_Phantom_4_regions_AO_7W_L5_100b.Rdata")
#load("./Data/results_Chicago_1024_9_AO_200b_36L.Rdata")#OK
#load("./Data/results_Mexico_512_7_AO_300b.Rdata")#OK
#load("./Data/results_panama_512_AO_9W_L5_100b.Rdata")
#save(difference_values, x, file = "./Data/results_Panama_512_AO_7x7_L5_200b.Rdata")
#load("./Data/results_Panama_512_AO_7x7_L5_200b.Rdata")
#load("./Data/results_Lake_512_9W_AO_100b_36L.Rdata")#OK
#load("./Data/results_Phantom_4_z1_200b.Rdata")#OK
#save(difference_values, Z, file = "./Data/results_Phantom_4_renyi_10.Rdata")
#save(difference_values, Z, file = "./Data/results_Phantom_4_renyi_B100.Rdata")
#load("./Data/results_Phantom_4_renyi_B100.Rdata")
#save(difference_values, x, file = "./Data/results_rotterdam_renyi_w7_10_07.Rdata")
#load("./Data/results_rotterdam_renyi_w7_10_07.Rdata")
#save(p_values_beta1, Z, file = "./Data/p_values_beta1.Rdata")
#save(p_values_beta1, Z, file = "./Data/p_values_beta1_v1.Rdata")
#load("./Data/p_values_beta1_v1.Rdata")
#save(difference_values, x, file = "./Data/results_rotterdam_renyi_w7_b200_0_65.Rdata")
#load("./Data/results_rotterdam_renyi_w7_b200_0_65.Rdata")
#load("./Data/results_Rotterdam_AO_L1_7_2.Rdata")
#save(difference_values, Z, file = "./Data/results_Phantom_4_renyi_B100_w7_075_L5.Rdata")#Ok but w9
#load("./Data/results_Phantom_4_renyi_B100_w7_075_L5.Rdata")
#save(difference_values, x, file = "./Data/results_panama_512_renyi_w7_b100_075_L5_w9.Rdata")
#load("./Data/results_panama_512_renyi_w7_b100_075_L5_w9.Rdata")
#save(difference_values, x, file = "./Data/results_illinois_1024_renyi_b100_085_L36_w9.Rdata")
#load("./Data/results_illinois_1024_renyi_b100_085_L36_w9.Rdata")
#save(difference_values, x, file = "./Data/results_sanFrancisco_650_renyi_b200_075_L5_w11.Rdata")
#load("./Data/results_sanFrancisco_650_renyi_b200_075_L5_w11.Rdata")
#save(difference_values, x, file = "./Data/results_illinois_1024_renyi_b200_075_L36_w11.Rdata")
#load("./Data/results_illinois_1024_renyi_b200_075_L36_w11.Rdata")
#save(difference_values, y, file = "./Data/results_rotterdam_renyi_w9_099.Rdata")
#load("./Data/results_rotterdam_renyi_w9_099.Rdata")
#load("./Data/results_rotterdam_renyi_w11_099.Rdata")
#load("./Data/results_rotterdam_new2_renyi_w11_099.Rdata")
#load("./Data/results_rotterdam_new2_Ebrahimi_w11.Rdata")
#save(difference_values, y, file = "./Data/results_rotterdam_origin_Ebrahimi_w11.Rdata")
load("./Data/results_rotterdam_origin_Ebrahimi_w11.Rdata")
#load("./Data/results_rotterdam_origin_Ebrahimi_w7.Rdata")
calculate_p_values_matrix <- function(data_matrix, mu, sigma) {
  rows <- nrow(data_matrix)
  cols <- ncol(data_matrix)
  
  p_values_matrix <- matrix(NA, nrow = rows, ncol = cols)
  
  for (i in 1:rows) {
    for (j in 1:cols) {
      test_difference <- data_matrix[i, j]
      
      epsilon <-(test_difference - 0) / (sigma)
      p_value <- 2 * pnorm(-abs(epsilon))#  2 * (1 - pnorm(abs(epsilon)))
      #2*pnorm(abs(epsilon))-1
      
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
#hist(p_values_matrix)
#plot(imagematrix(p_values_matrix),significance_level = 0.05)
#plot(imagematrix(p_values_matrix ))
plot(imagematrix(p_values_matrix >0.05))
#plot(imagematrix(equalize(difference_values)))
#plot(imagematrix(equalize(x)))


#      imagematrixPNG(imagematrix(equalize(y)), name = "rotterdam_new2_1024.png")
#   imagematrixPNG(imagematrix(equalize(difference_values)), name = "rotterdam_ebrahimi_w11.png")
# # # # # # # # # # # 
#   imagematrixPNG(imagematrix(p_values_matrix), name="H_pvalue_rotterdam_new2_1024_L1_w11_ebrahimi.png")
#   imagematrixPNG(imagematrix(p_values_matrix>0.05), name="H_005_rotterdam_origin_1024_L1_w11_ebrahimi.png")
