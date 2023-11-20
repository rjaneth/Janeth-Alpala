# Function for Generating Samples from the GI0 Distribution

# Input:
#   - L: Looks.
#   - alpha: Roughness parameter.
#   - mu: The mean parameter 
#   - n: Sample size
#
# Output:
#   - Random samples following the GI0 distribution.



gi0_sample <- function(mu, alpha, L,  n) {
  
  X_samples <- rinvgamma(n, shape = -alpha, rate =mu*(-alpha-1))
  
  Y_samples <- rgamma(n, shape = L, rate = L)
  Z_samples <- X_samples * Y_samples
  return(Z_samples)
}