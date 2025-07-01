#Comparacion de estimadores

# --- estimador rápido (sin loops explícitos) ---------------------------
renyi_entropy_spacing_1 <- function(x, lambda, m = floor(sqrt(length(x)) + 0.5))
{
  n <- length(x)
  x <- sort(x, method = "quick")
  i <- seq_len(n)
  
  li <- pmax(i - m, 1L)
  ri <- pmin(i + m, n)
  diff <- x[ri] - x[li]
  
  ci <- ifelse(i <= m, (m + i - 1)/m,
               ifelse(i >= n - m + 1, (n + m - i)/m, 2))
  
  r <- n * diff / (ci * m)
  r[r <= .Machine$double.eps] <- NA
  s <- mean(r^(1 - lambda), na.rm = TRUE)
  
  (1 / (1 - lambda)) * log(s)
}

# --- bootstrap acelerado ----------------------------------------------
bootstrap_renyi_entropy <- function(x, B = 200L, lambda,
                                    m = floor(sqrt(length(x)) + 0.5),
                                    parallel = FALSE)
{
  n <- length(x)
  
  # estimador en los datos originales
  t0 <- renyi_entropy_spacing(x, lambda, m)
  
  # matriz de índices bootstrap (n × B)
  idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n)
  
  # función auxiliar (usa la misma m para todas las réplicas)
  f_boot <- function(col) renyi_entropy_spacing(x[col], lambda, m)
  
  # vector de estimaciones bootstrap
  if (parallel) {
    boot_vals <- future.apply::future_apply(idx, 2, f_boot)
  } else {
    boot_vals <- apply(idx, 2, f_boot)
  }
  
  # sesgo estimado  = mean(bootstrap) - t0
  bias_hat <- mean(boot_vals, na.rm = TRUE) - t0
  
  # versión “bias-corrected”   2·t0 − mean(bootstrap)
  2 * t0 - mean(boot_vals, na.rm = TRUE)
}


renyi_entropy_estimator_v1 <- function(data, a) {
  n <- length(data)
  m <- round(sqrt(n) + 0.5)  # Espaciado m
  data_sorted <- sort(data)   # Ordenar los datos
  
  # Calcular c_i para cada valor de i usando vectores
  c_i <- numeric(n)
  c_i[1:m] <- 1 + (0:(m-1)) / m
  c_i[(m+1):(n-m)] <- 2
  c_i[(n-m+1):n] <- 1 + (n - (n-m+1):n) / m
  
  # Manejar los límites y calcular diff_term para cada i
  diff_term <- numeric(n)
  for (i in 1:n) {
    if (i <= m) {
      diff_term[i] <- data_sorted[i + m] - data_sorted[1]
    } else if (i >= n - m + 1) {
      diff_term[i] <- data_sorted[n] - data_sorted[i - m]
    } else {
      diff_term[i] <- data_sorted[i + m] - data_sorted[i - m]
    }
  }
  
  # Calcular la suma para el estimador
  sum_term <- sum(((n / (c_i * m)) * diff_term)^(1 - a))
  
  # Calcular la entropía de Rényi usando la fórmula general
  renyi_entropy <- (1 / (1 - a)) * log(sum_term / n)
  
  return(renyi_entropy)
}

bootstrap_renyi_entropy_estimator_v1 <- function(x, B, a) {
  
  n <- length(x)
  v.Bootstrap <- numeric(B)
  
  for (b in 1:B) {
    # Muestra bootstrap con reemplazo
    bx <- sample(x, size = n, replace = TRUE)
    
    # Verificar si todos los valores no son iguales para evitar un resultado no informativo
    if (length(unique(bx)) > 1) {
      entropy_result <- renyi_entropy_estimator_v1(bx, a)
      if (is.finite(entropy_result)) {
        v.Bootstrap[b] <- entropy_result
      } else {
        v.Bootstrap[b] <- NA
      }
    } else {
      v.Bootstrap[b] <- NA  # Asignar NA si todos los valores son iguales
    }
  }
  
  # Calcular la entropía de Rényi para la muestra original
  t <- renyi_entropy_estimator_v1(x, a)
  
  # Calcular el promedio de los valores bootstrap válidos
  estimated_mean <- mean(v.Bootstrap, na.rm = TRUE)
  
  # Estimador mejorado
  return(2 * t - estimated_mean)
}

# library(microbenchmark)
# set.seed(1)
# x <- rgamma(1e4, shape = 2)
# lambda <- 1.5
# B <- 100
# 
# microbenchmark(
#   bootstrap_slow = bootstrap_renyi_entropy_estimator_v1(x, B, lambda), # tu versión
#   bootstrap_fast = bootstrap_renyi_entropy        (x, B, lambda),      # nueva versión
#   times = 5
# )

library(microbenchmark)
set.seed(1)

n <- 50000      # muestra más grande
B <- 200      # más repeticiones
lambda <- 1.5
x <- gi0_sample(3, -3.5, 4, n)

microbenchmark(
  slow = bootstrap_renyi_entropy_estimator_v1(x, B, lambda),
  fast = bootstrap_renyi_entropy(x, B, lambda),
  times = 3
)

# # Necesarios
# library(invgamma)
# library(microbenchmark)
# 
# # Si vas a usar paralelización más adelante
# # install.packages("future.apply")
# gi0_sample <- function(mu, alpha, L, n) {
#   if (alpha >= -1) stop("alpha debe ser < -1")
#   if (mu <= 0 || L <= 0) stop("mu y L deben ser > 0")
#   
#   X <- rinvgamma(n, shape = -alpha, rate = mu * (-alpha - 1))
#   Y <- rgamma(n, shape = L, rate = L)
#   Z <- X * Y
#   return(Z)
# }
# set.seed(123)
# n <- 10000
# B <- 1000
# lambda <- 1.5
# mu <- 3
# alpha <- -3.5
# L <- 4
# 
# # Simulación de muestra GI0
# x <- gi0_sample(mu, alpha, L, n)
# 
# # Comparar resultados
# res_slow <- bootstrap_renyi_entropy_estimator_v1(x, B, lambda)
# res_fast <- bootstrap_renyi_entropy(x, B, lambda)
# 
# cat("Resultado original (v1): ", res_slow, "\n")
# cat("Resultado optimizado   : ", res_fast, "\n")
# cat("Diferencia absoluta    : ", abs(res_slow - res_fast), "\n")
