#571 
# L <- 5
# s <- 49
# R <- 100000
# mu <- 1
# 
# TestStat1 <- numeric(R)
# 
# for (r in 1:R) {
#   z <- gamma_sar_sample(L, mu, s)
#   Q2 <- median(z)
#   MnAD <- mean(abs(z - Q2))
#   TestStat1[r] <- MnAD/Q2
# }
# 
# mean(TestStat1) # Valor esperado bajo homogeneidad


L <- 5
s <- 49
R <- 10000
mu <- 1

TestStat1 <- numeric(R)

for (r in 1:R) {
  z <- gamma_sar_sample(L, mu, s) # Genera datos homogéneos
  Q2 <- median(z) # Mediana muestral
  MnAD <- mean(abs(z - Q2)) # MnAD respecto a la mediana
  TestStat1[r] <- MnAD / Q2 # Estadístico inicial
}

mean_expected <- mean(TestStat1) # Valor esperado bajo homogeneidad
print(mean_expected)
TestStat1_adjusted <- numeric(R)

for (r in 1:R) {
  z <- gamma_sar_sample(L, mu, s) # Genera datos homogéneos
  Q2 <- median(z) # Mediana muestral
  MnAD <- mean(abs(z - Q2)) # MnAD respecto a la mediana
  TestStat1_adjusted[r] <- (MnAD / Q2) - mean_expected # Ajustar
}

# Validación: El promedio del estadístico ajustado debe ser cercano a cero
print(mean(TestStat1_adjusted)) # Debería ser cercano a 0

