gi0_density <- function(x, mu, alpha, L) {
  
  X_samples <- rinvgamma(n, shape = -alpha, rate =mu*(-alpha-1))
  
  Y_samples <- rgamma(n, shape = L, rate = L)
  Z_samples <- X_samples * Y_samples
  return(Z_samples)
}
