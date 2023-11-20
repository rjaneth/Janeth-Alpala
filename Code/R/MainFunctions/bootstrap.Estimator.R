#
source("../MainFunctions/van_es_estimator.r")
source("../MainFunctions/correa_estimator.r")
source("../MainFunctions/ebrahimi_estimator.r")
source("../MainFunctions/noughabi_arghami_estimator.r")
source("../MainFunctions/vasicek_estimator.r")
source("../MainFunctions/al_omari_1_estimator.r")
source("../MainFunctions/al_omari_2_estimator.r")

# bootstrap.Estimator <- function(x, B){
#   sample_size <- length(x)
#   t <- van_es_estimator(x)
#   
#   v.Bootstrap <- rep(0, B) 
#   for(b in 1:B) {
#     bx <- sample(x, replace = TRUE, size = sample_size)  # Toma muestras del mismo tamaño que la muestra original
#     
#     # Añade una condición para verificar si todos los valores son iguales en la muestra Bootstrap
#     if (!all(bx == bx[1])) {
#       est_bx <- van_es_estimator(bx)
#       
#       # Añade una condición para descartar estimaciones donde la entropía es cero
#       if (est_bx == 0) {
#         v.Bootstrap[b] <- NA
#       } else {
#         v.Bootstrap[b] <- est_bx
#       }
#     } else {
#       v.Bootstrap[b] <- NA  # Descarta la muestra Bootstrap con todos los valores iguales
#     }
#   }
#   
#   bootstrap_mean <- mean(v.Bootstrap, na.rm = TRUE)
#   
#   improved_estimator <- 2 * t - bootstrap_mean
#   
#   return(improved_estimator)
# }



# bootstrap.Estimator <- function(x, B){
#   sample_size <- length(x)
#   t <- van_es_estimator(x)
#   
#   # Define B en función del tamaño de la muestra
#   #B <- sample_size * 2  # Puedes ajustar esto según tus necesidades
#   
#   v.Bootstrap <- rep(0, B) 
#   for(b in 1:B) {
#     bx <- sample(x, replace = TRUE, size = sample_size)  # Toma muestras del mismo tamaño que la muestra original
#     est_bx <- van_es_estimator(bx)
#     
#     # Puedes agregar aquí una condición para descartar si la entropía es cero
#     if (est_bx == 0) {
#       v.Bootstrap[b] <- NA
#     } else {
#       v.Bootstrap[b] <- est_bx
#     }
#   }
#   
#   bootstrap_mean <- mean(v.Bootstrap, na.rm = TRUE)
#   
#   improved_estimator <- 2 * t - bootstrap_mean
#   
#   return(improved_estimator)
# }

bootstrap.Estimator <- function(x, B){
  sample_size <- length(x)
  t <- van_es_estimator(x)

  v.Bootstrap <- rep(0, B)
  for(b in 1:B) {
    bx <- sample(x, replace = TRUE, size = sample_size)
    est_bx <- van_es_estimator(bx)
    v.Bootstrap[b] <- est_bx


  #   if (est_bx == 0) {
  #     v.Bootstrap[b] <- NA
  #   }
 }


  bootstrap_mean <- mean(v.Bootstrap, na.rm = TRUE)


  improved_estimator <- 2 * t - bootstrap_mean

  return(improved_estimator)
}
