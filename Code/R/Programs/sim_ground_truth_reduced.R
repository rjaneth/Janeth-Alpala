#559
# Cargar librerías necesarias
library(png)
source("../imagematrix.R")

# Leer la imagen base Phantom1 (binaria en blanco y negro)
Phantom1 <- readPNG("../../../Figures/PNG/Phantom1.png")[,,1]
dim(Phantom1)  # Confirmar dimensiones (500x500)

# Crear matriz binaria del ground truth (ajustada para tamaño 500x500)
ground_truth <- matrix(0, nrow = 500, ncol = 500)

# Definir regiones según las máscaras de la imagen Phantom1
# Zona superior izquierda
p.up.le <- Phantom1[1:250, 1:250]
ground_truth[1:250, 1:250][p.up.le == 1] <- 1  # Región heterogénea (α = -1.5)

# Zona superior derecha
p.up.ri <- Phantom1[1:250, 251:500]
ground_truth[1:250, 251:500][p.up.ri == 1] <- 1  # Región heterogénea (α = -2)

# Zona inferior izquierda
p.bo.le <- Phantom1[251:500, 1:250]
ground_truth[251:500, 1:250][p.bo.le == 1] <- 1  # Región heterogénea (α = -2.5)

# Zona inferior derecha
p.bo.ri <- Phantom1[251:500, 251:500]
ground_truth[251:500, 251:500][p.bo.ri == 1] <- 1  # Región heterogénea (α = -3)

# Reducir la matriz binaria para que coincida con las dimensiones de las ventanas deslizantes (494x494)
ground_truth_reduced <- ground_truth[4:497, 4:497]

# Guardar las matrices binaria completas y reducidas
#save(ground_truth, file = "./Data/ground_truth_binary_full.Rdata")
#save(ground_truth_reduced, file = "./Data/ground_truth_binary_reduced.Rdata")

# Visualizar las matrices
#par(mfrow = c(1, 2))
#plot(imagematrix(ground_truth))
plot(imagematrix(ground_truth_reduced))
