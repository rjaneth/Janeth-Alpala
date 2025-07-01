
source("./Code/entropy_gamma_sar.R")
source("./Code/entropy_gI0.R")
source("./Code/gamma_sar_sample.R")
source("./Code/gi0_sample.R")
source("./Code/al_omari_1_estimator.R")
source("./Code/bootstrap_al_omari_1_estimator.R")
source("./Code/renyi_entropy_estimator_v1.R")
source("./Code/bootstrap_renyi_entropy_estimator_v1.R")
source("./Code/entropy_renyi_gamma_sar.R")
source("./Code/functions_sample_bias_mse.R")
source("./Code/functions_sample_bias_mse_1.R")

source("./Code/tsallis_entropy_gamma_sar.R")
source("./Code/tsallis_estimator.R")
source("./Code/bootstrap_tsallis_entropy.R")
source("./Code/tsallis_estimator_optimized.R")
source("./Code/bootstrap_tsallis_entropy_optimized.R")

source("./Code/renyi_estimator.R")
source("./Code/bootstrap_renyi_estimator.R")
source("./Code/renyi_entropy_optimized.R")
source("./Code/bootstrap_renyi_estimator_Op.R")

source("./Code/shannon_alomari_estimator.R")
source("./Code/bootstrap_shannon_alomari.R")
source("./Code/shannon_entropy_gamma_sar.R")

#------------------------------------------------------------
# 2.  Bootstrap bias-corrected para Tsallis
#------------------------------------------------------------
bootstrap_tsallis_entropy <- function(x, B = 200L, alpha,
                                      m = floor(sqrt(length(x)) + 0.5),
                                      parallel = FALSE)
{
  n  <- length(x)
  
  # estimación en la muestra original
  t0 <- tsallis_estimator(x, alpha, m)
  
  # matriz de índices bootstrap (n × B)
  idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n)
  
  # función auxiliar (m mismo p/ todas las réplicas)
  f_boot <- function(col) tsallis_estimator(x[col], alpha, m)
  
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