tsallis_gammasar_direct <- function(mu, L, lambda) {
  # Validación de parámetros
  if (mu <= 0 || L <= 0 || lambda <= 0) {
    stop("Parámetros deben ser positivos: mu > 0, L > 0, lambda > 0")
  }
  
  # Calcular términos individuales
  term1 <- mu^(1 - lambda)
  term2 <- L^(lambda - 1)
  term3 <- lambda^(-lambda*(L-1) - 1)
  term4 <- gamma(lambda*(L-1) + 1)
  term5 <- gamma(L)^(-lambda)
  
  # Calcular el término completo dentro del corchete
  bracket_term <- term1 * term2 * term3 * term4 * term5
  
  # Entropía de Tsallis
  (1 - bracket_term) / (lambda - 1)
}

#expresiones compactas teoricas

tsallis_gammasar <- function(mu, L, lambda)
  (1 - mu^(1-lambda) *
     L^(lambda-1) *
     lambda^(-lambda*(L-1)-1) *
     gamma(lambda*(L-1)+1) *
     gamma(L)^(-lambda)
  ) / (lambda - 1)

tsallis_gi0 <- function(alpha, mu, L, lambda)
  (1 - mu^(1-lambda) *                     # μ^{1-λ}
     L^(lambda-1) *                      # L^{λ-1}
     (-alpha - 1)^(1-lambda) *           # (−α−1)^{1-λ}
     (gamma(L - alpha) / gamma(-alpha))^lambda *  # [Γ(L-α)/Γ(−α)]^λ
     gamma(lambda*(L-1) + 1) *           # Γ(λ(L−1)+1)
     gamma(L)^(-lambda) *                # Γ(L)^{−λ}
     gamma(lambda*(1 - alpha) - 1) /     # Γ(λ(1−α)−1)
     gamma(lambda*(L - alpha))           # Γ(λ(L−α))
  ) / (lambda - 1)


# --- Término base Γ_SAR(L,μ) en forma logarítmica (una sola línea) ----
tsallis_gammasar_log <- function(mu, L, lambda)
  (1 - exp((1 - lambda)*log(mu) +
             (lambda - 1)*log(L) +
             lgamma(lambda*(L - 1) + 1) -
             lambda*lgamma(L) -
             (lambda*(L - 1) + 1)*log(lambda))) /
  (lambda - 1)

tsallis_entropy_gamma_sar <- function(mu, L, lambda) {
  (1 - exp((1 - lambda)*log(mu) +
             (lambda - 1)*log(L) +
             lgamma(lambda*(L - 1) + 1) -
             lambda*lgamma(L) -
             (lambda*(L - 1) + 1)*log(lambda))) / (lambda - 1)
}

# --- Entropía total T_λ(GI0) en forma log-exp (una sola línea) --------
tsallis_gi0_log <- function(alpha, mu, L, lambda)
  (1 - exp((1 - lambda)*log(mu) +
             (lambda - 1)*log(L) +
             (1 - lambda)*log(-alpha - 1) +                     # nuevo factor textura
             lgamma(lambda*(L - 1) + 1) -
             lambda*lgamma(L) +
             lambda*(lgamma(L - alpha) - lgamma(-alpha)) +       # cociente Γ(L-α)/Γ(−α)
             lgamma(lambda*(1 - alpha) - 1) -
             lgamma(lambda*(L - alpha)))) /
  (lambda - 1)
