# functions_bias_mse.R

generate_samples <- function(sample_size, replication, mu, L) {
  samples <- vector("list", replication)
  for (r in 1:replication) {
    samples[[r]] <- gamma_sar_sample(L, mu, sample_size)
  }
  return(samples)
}

calculate_bias_mse <- function(sample_sizes, R, mu, L) {
  true_entropy <- entropy_gamma_sar(L, mu)
  
  # Resto del código...
  
  return(output)
}
