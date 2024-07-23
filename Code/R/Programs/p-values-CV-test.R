
#rm(list = ls())
#save( cv_values, x, file = "./Data/CV_results_data_Illinois_crops_1024.Rdata")
#load("../Programs/Data/CV_results_data_Illinois_crops_1024.Rdata")

#load("../Programs/Data/results_data_simulated_Phantom_7.Rdata")

#load("../Programs/Data/CV_results_data_Mexixo_512.Rdata")
#load("../Programs/Data/CV_results_data_lake_512.Rdata")
#load("../Programs/Data/MnAD_results_data_lake_512.Rdata")
#load("../Programs/Data/CV_results_data_Phantom_4_z1.Rdata")
load("../Programs/Data/CV_results_data_Rotterdam_1024.Rdata")

# Cargar los datos y parámetros de ajuste de la distribución lognormal
#Illinois_crops_1024
# meanlog <- -1.3877655
# sdlog <-  0.5221934


#MEXICO
# 
# meanlog <- -1.0042160
# sdlog <-  0.4241594

#lake

 # meanlog <- -1.1116401
 # sdlog <-  0.5720475

# #phantom z1
# meanlog <- 0.0164172
# sdlog <-   0.5601495

#Rottember
meanlog <- -0.1182653
sdlog <-   0.3187516

# Crear una matriz para almacenar los p-valores
p_values_matrix <- matrix(NA, nrow = nrow(cv_values), ncol = ncol(cv_values))

# Calcular los p-valores para cada dato en la matriz de coeficientes de variación
for (i in 1:nrow(cv_values)) {
  for (j in 1:ncol(cv_values)) {
    cv_value <- cv_values[i, j]
    
    # Calcular el p-valor utilizando la función plnorm()
    p_values_matrix[i, j] <- 1-plnorm(cv_value, meanlog = meanlog, sdlog = sdlog)
  }
}


#save(p_values_matrix, file = "./Data/results_pvalue_Flevoland_cv_300_5.Rdata")




#source("../imagematrix.R")
source("../imagematrix.R")
#hist(p_values_matrix)
#par(mfrow=c(1,2))

#plot(imagematrix(equalize(cv_values)))#z.up.le
plot(imagematrix(p_values_matrix ))
#plot(imagematrix(p_values_matrix >0.1))
plot(imagematrix(p_values_matrix >0.05))

# # #imagematrixPNG(imagematrix(equalize(x)), name = "Phantom_AO_7W_L5_100b.png")
 imagematrixPNG(imagematrix(equalize(cv_values)), name = "cv_Rotterdam.png")
# # 
 imagematrixPNG(imagematrix(p_values_matrix), name="cv_pvalues_Rotterdam.png")
  imagematrixPNG(imagematrix(p_values_matrix>0.05), name="cv_005_pvalues_Rotterdam.png")
