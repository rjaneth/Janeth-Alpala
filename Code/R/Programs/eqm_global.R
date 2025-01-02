
library(fields)


load("./Data/p_values_renyi_4z_renyi_B200_w7_09_L5.Rdata")


load("./Data/p_values_shannon_4z_renyi_B200_w7_09_L5.Rdata")

#  mapas binarios con umbral de significancia (alpha = 0.05)
alpha <- 0.05  

binary_map_shannon <- ifelse(p_values_shannon <= alpha, 1, 0)
#save(binary_map_shannon, file = "./Data/binary_map_shannon.Rdata")
binary_map_renyi <- ifelse(p_values_renyi <= alpha, 1, 0)
#save(binary_map_renyi, file = "./Data/binary_map_renyi.Rdata")

# EQM Local (pixel por pixel)
eqm_local_shannon <- (binary_map_shannon - ground_truth_reduced)^2
eqm_local_renyi <- (binary_map_renyi - ground_truth_reduced)^2

# plot
par(mfrow = c(1, 2))


# Shannon
image.plot(eqm_local_shannon, main = "Errores Locales - Shannon", 
           xlab = "Columnas", ylab = "Filas", col = terrain.colors(100))

#  Rényi
image.plot(eqm_local_renyi, main = "Errores Locales - Rényi", 
           xlab = "Columnas", ylab = "Filas", col = terrain.colors(100))

#   EQM Global a partir del EQM Local
eqm_global_shannon <- mean(eqm_local_shannon, na.rm = TRUE)
eqm_global_renyi <- mean(eqm_local_renyi, na.rm = TRUE)

# 
cat("EQM Global - Shannon:", eqm_global_shannon, "\n")
cat("EQM Global - Rényi:", eqm_global_renyi, "\n")