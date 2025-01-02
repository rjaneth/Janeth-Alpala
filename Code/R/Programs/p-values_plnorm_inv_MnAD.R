library(ggplot2)

library(reshape2)
#library(plotly)
library(knitr)
library(pandoc)
library(gridExtra)
library(tidyr)
library(gtools)
library(stats4)
library(rmutil)
library(scales)
library(tidyr)
library(gtools)
library(stats4)
library(rmutil)
library(invgamma)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ggsci)
library(wesanderson)
library(ggpubr)
library(patchwork)
#options(kableExtra.latex.load_packages = FALSE)
library(devtools)
#devtools::install_github("haozhu233/kableExtra")
#devtools::install_github("haozhu233/kableExtra")
library(kableExtra)
library(ggthemes)
library(geomtextpath)
library(latex2exp)
library(univariateML)
library(viridis)
#install.packages('univariateML')
#univariateML::univariateML_models

library(fitdistrplus)
library(logspline)
theme_set(theme_bw()  +
            theme(text=element_text(family="serif"),
                  legend.position = "top")# Gridtop , right , bottom , or left#, panel.grid = element_blank()
)
# Helpful for latex tables
# library(xtable)
# options(xtable.caption.placement='top',
#         xtable.table.placement='!t',
#         xtable.include.rownames=F,
#         xtable.comment=F)
# #library(rgl)



source("../MainFunctions/gamma_sar_sample.r")
source("../MainFunctions/entropy_gamma_sar.r")
source("../MainFunctions/entropy_gI0.r")

source("../MainFunctions/van_es_estimator.r")
source("../MainFunctions/correa_estimator.r")
source("../MainFunctions/ebrahimi_estimator.r")
source("../MainFunctions/noughabi_arghami_estimator.r")
source("../MainFunctions/vasicek_estimator.r")
source("../MainFunctions/al_omari_1_estimator.r")
source("../MainFunctions/al_omari_2_estimator.r")

source("../MainFunctions/bootstrap_van_es_estimator.r")
source("../MainFunctions/bootstrap_correa_estimator.r")
source("../MainFunctions/bootstrap_ebrahimi_estimator.r")
source("../MainFunctions/bootstrap_noughabi_arghami_estimator.r")
source("../MainFunctions/bootstrap_vasicek_estimator.r")
source("../MainFunctions/bootstrap_al_omari_1_estimator.r")
source("../MainFunctions/bootstrap_al_omari_2_estimator.r")
#The next function contains the functions: generate_samples, calculate_bias_mse, generate_plot
source("../Programs/functions_sample_bias_mse.R")# 
source("../imagematrix.R")

#save(cd_values_mnad, x, file = "./Data/results_data_Frankfurt_512_mnad_7.Rdata")
load("../Programs/Data/results_data_Frankfurt_512_mnad_7.Rdata")
#load("../Programs/Data/results_data_simulated_z_up_le_cv_250_7.Rdata")
# Definir la función de distribución acumulativa inversa de la distribución lognormal
plnorm_inv <- function(x, meanlog, sdlog) {
  pnorm(log(x), mean = meanlog, sd = sdlog)
}

# Cargar los datos y parámetros de ajuste de la distribución lognormal
meanlog <- -0.1931404
sdlog <- 0.4293747

# Crear una matriz para almacenar los p-valores
p_values_matrix <- matrix(NA, nrow = nrow(cv_values), ncol = ncol(cv_values))

# Calcular los p-valores para cada dato en la matriz de coeficientes de variación
for (i in 1:nrow(cv_values)) {
  for (j in 1:ncol(cv_values)) {
    cv_value <- cv_values[i, j]
    
    # Calcular el p-valor utilizando la función de distribución acumulativa inversa de la distribución lognormal
    p_values_matrix[i, j] <- plnorm_inv(cv_value, meanlog, sdlog)
  }
}

# Guardar la matriz de p-valores
#save(p_values_matrix, file = "./Data/results_pvalue_Flevoland_cv_300_5.Rdata")


#source("../imagematrix.R")
#source("../imagematrix.R")
#hist(p_values_matrix)
#par(mfrow=c(1,2))
#plot(imagematrix(p_values_matrix ))
plot(imagematrix(equalize(z.up.le)))
plot(imagematrix(equalize(cv_values)))#z.up.le
plot(imagematrix(p_values_matrix ))
plot(imagematrix(p_values_matrix <0.1))
plot(imagematrix(p_values_matrix >0.1))

imagematrixPNG(imagematrix(equalize(x)), name = "simulated7x7.png")
imagematrixPNG(imagematrix(p_values_matrix>0.01), name="sim_pvalue7x7.png")