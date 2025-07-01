# no se si funciona
renyi_entropy_spacing <- function(data, lambda, m = NULL)
{
  # Validaciones básicas
  if (lambda == 1) stop("Para λ = 1 usa el estimador de Shannon")
  n <- length(data)
  if (n < 10) warning("Muestra muy pequeña; resultados inestables")
  if (is.null(m)) m <- max(1L, round(sqrt(n) + 0.5))
  
  x  <- sort(data, method = "quick")
  i  <- seq_len(n)
  
  # índices de los extremos del m-spacing
  li <- pmax(i - m, 1L)      # left  index
  ri <- pmin(i + m, n)       # right index
  
  # diferencias de forma vectorizada
  diff_term <- x[ri] - x[li]
  
  # coeficiente c_i con la misma lógica de bordes
  c_i <- ifelse(i <= m,
                (m + i - 1)/m,
                ifelse(i >= n - m + 1,
                       (n + m - i)/m,
                       2))
  
  # razón r_i = n/(c_i m) * diff
  r <- n * diff_term / (c_i * m)
  r[r <= .Machine$double.eps] <- NA               # evita log(0)
  
  # promedio de r_i^{1-λ}
  s <- mean(r^(1 - lambda), na.rm = TRUE)
  
  (1 / (1 - lambda)) * log(s)
}
