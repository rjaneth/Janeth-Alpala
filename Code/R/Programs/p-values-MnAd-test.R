rm(list = ls())

#load("../Programs/Data/results_Phantom_mnad_7.Rdata")
#save(cd_values_mnad, Z, file = "./Data/results_Phantom_mnad_7.Rdata")

#load("../Programs/Data/MnAD_results_data_Illinois_crops_1024.Rdata")
#save( cd_values_mnad, x, file = "./Data/MnAD_results_data_mexico_512.Rdata")
#load("../Programs/Data/MnAD_results_data_mexico_512.Rdata")
#load("../Programs/Data/MnAD_results_data_lake_512.Rdata")
#load("../Programs/Data/results_Phantom_mnad_7_z1.Rdata")
load("../Programs/Data/results_Phantom_mnad_7_new.Rdata")

#Illinois_crops
# meanlog <- -1.5865643
# sdlog <- 0.5912683
#mexico
# meanlog <- -1.1788140
# sdlog <- 0.5000859

# #lake
# meanlog <- -1.2856986
# sdlog <- 0.6665039

# #phantom z1
#  meanlog <- -0.005317387
#  sdlog <- 0.783799230
  meanlog <- 0.004045258
  sdlog <- 0.848081854

# Crear una matriz para almacenar los p-valores
p_values_matrix <- matrix(NA, nrow = nrow(cd_values_mnad), ncol = ncol(cd_values_mnad))

# Calcular los p-valores para cada dato en la matriz de coeficientes de variación
for (i in 1:nrow(cd_values_mnad)) {
  for (j in 1:ncol(cd_values_mnad)) {
    cd_value_mnad <- cd_values_mnad[i, j]
    
    # Calcular el p-valor utilizando la función plnorm()
    p_values_matrix[i, j] <- 1-plnorm(cd_value_mnad, meanlog = meanlog, sdlog = sdlog)
  }
}

# Guardar la matriz de p-valores
#save(p_values_matrix, file = "./Data/results_pvalue_Flevoland_cv_300_5.Rdata")




#source("../imagematrix.R")
source("../imagematrix.R")
hist(cd_values_mnad)
#par(mfrow=c(1,2))
#plot(imagematrix(p_values_matrix ))
#plot(imagematrix(equalize(z.up.le)))
plot(imagematrix(equalize(cd_values_mnad)))#z.up.le
plot(imagematrix(p_values_matrix ))
#plot(imagematrix(p_values_matrix <0.1))
plot(imagematrix(p_values_matrix >0.05))

# imagematrixPNG(imagematrix(equalize(Z)), name = "mnad_Phantom_mnad_7_z1.png")
# 
# imagematrixPNG(imagematrix(equalize(cd_values_mnad)), name = "Phantom_mnad_7_z1.png")
# # # # 
# imagematrixPNG(imagematrix(p_values_matrix), name="mnad_p_values_Phantom_mnad_7_z1.png")
# imagematrixPNG(imagematrix(p_values_matrix>0.05), name="mnad_005_Phantom_7_z1.png")
