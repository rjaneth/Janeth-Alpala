#------------------------------------------------------------
# 1.  Estimador base (spacing) de Tsallis  –  vectorizado
#------------------------------------------------------------
# comparacion 2
tsallis_entropy_spacing2 <- function(x, alpha,
                                    m = floor(sqrt(length(x))+ 0.5 ))
{
  if (alpha <= 0) stop("alpha debe ser > 0 (Tsallis)")
  n <- length(x)
  if (n < 10) warning("Muestra muy pequeña; resultados inestables")
  
  xs <- sort(x, method = "quick")
  i  <- seq_len(n)
  
  li <- pmax(i - m, 1L)         # índices i-m (cortados a 1)
  ri <- pmin(i + m, n)          # índices i+m (cortados a n)
  
  diff <- xs[ri] - xs[li]
  
  ci <- ifelse(i <= m,             (m + i - 1)/m,
               ifelse(i >= n - m + 1,     (n + m - i)/m, 2))
  
  r <- n * diff / (ci * m)
  r[r <= .Machine$double.eps] <- NA
  
  s <- sum(r^(1 - alpha), na.rm = TRUE)  # ∑ r_i^{1-α}
  
  (1 / (alpha - 1)) * (1 - s / n)
}

#------------------------------------------------------------
# 2.  Bootstrap bias-corrected para Tsallis
#------------------------------------------------------------
bootstrap_tsallis_entropy <- function(x, B = 200L, alpha,
                                      m = floor(sqrt(length(x)) + 0.5),
                                      parallel = FALSE)
{
  n  <- length(x)
  
  # estimación en la muestra original
  t0 <- tsallis_entropy_spacing(x, alpha, m)
  
  # matriz de índices bootstrap (n × B)
  idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n)
  
  # función auxiliar (m mismo p/ todas las réplicas)
  f_boot <- function(col) tsallis_entropy_spacing(x[col], alpha, m)
  
  # estimaciones bootstrap
  if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE))
      stop("Instala el paquete 'future.apply' o usa parallel = FALSE")
    boot_vals <- future.apply::future_apply(idx, 2, f_boot)
  } else {
    boot_vals <- apply(idx, 2, f_boot)
  }
  
  # estimador bias-corregido:  2·t0 − mean(bootstrap)
  2 * t0 - mean(boot_vals, na.rm = TRUE)
}
tsallis_entropy_estimator1 <- function(data, alpha, m = NULL) {
  # Validación de parámetros
  if (alpha <= 0) stop("alpha debe ser mayor que 0")
  if (length(data) < 10) warning("Muestra muy pequeña, resultados pueden ser inestables")
  
  n <- length(data)
  if (is.null(m)) m <- max(1, round(sqrt(n)))
  data_sorted <- sort(data)
  
  # Cálculo vectorizado de diferencias
  left_idx <- pmax(1:n - m, 1)
  right_idx <- pmin(1:n + m, n)
  diff_term <- data_sorted[right_idx] - data_sorted[left_idx]
  
  # Factores de corrección de bordes
  c_i <- ifelse(
    1:n <= m, (m + 1:n - 1)/m,
    ifelse(
      1:n >= (n - m + 1), (n + m - 1:n)/m, 
      2
    )
  )
  
  # Cálculo de ratio con protección contra ceros
  ratio <- (n * diff_term) / (c_i * m)
  ratio[ratio <= .Machine$double.eps^0.5] <- NA
  
  # Término exponencial
  exp_term <- ratio^(1 - alpha)
  sum_exp <- sum(exp_term, na.rm = TRUE)
  n_valid <- sum(!is.na(exp_term))
  
  if (n_valid == 0) return(NA_real_)
  
  (1 - sum_exp/n) / (alpha - 1)
}
#------------------------------------------------------------
# 3.  Ejemplo rápido con datos GI0
#------------------------------------------------------------
library(invgamma)

gi0_sample <- function(mu, alpha, L, n)
{
  if (alpha >= -1) stop("alpha debe ser < -1")
  X <- rinvgamma(n, shape = -alpha, rate = mu * (-alpha - 1))
  Y <- rgamma(n, shape = L, rate = L)
  X * Y
}

set.seed(123)
n <- 2000
B <- 100
alpha_tsallis <- 0.7
mu <- 3; alpha_gi0 <- -300; L <- 14

x <- gi0_sample(mu, alpha_gi0, L, n)

# estimación Tsallis bootstrap (sin y con paralelización)
nnn <- tsallis_entropy_estimator1(x,alpha_tsallis)
t_bias_corrected  <- tsallis_entropy_spacing (x,  alpha_tsallis)
orig <- bootstrap_tsallis_entropy(x, B, alpha_tsallis)
# library(future); plan(multisession)             # para usar varios núcleos
# t_bias_corrected <- bootstrap_tsallis_entropy(x, B, alpha_tsallis, parallel = TRUE)
nnn
orig
t_bias_corrected
