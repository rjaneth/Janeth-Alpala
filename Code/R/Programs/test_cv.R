#570 test cv
set.seed(1234567890, kind = "Mersenne-Twister")

# Parámetros de simulación
R <- 10000
mu <- 1
sample.size <- c(49, 81, 121)
L_values <- c(3, 5, 8, 18)

all_summary_stats <- list()
all_TestStatistics <- list()

# # Función para calcular S_CV
# S_CV <- function(z, L) {
#   n <- length(z)
#   mean_z <- mean(z)
#   sd_z <- sd(z)
#   return((mean_z / sd_z) - (1 / sqrt(L)))
# }

# Simulación para cada valor de L y sample.size
for (L in L_values) {
  TestStatistics1 <- list()
  summary_stats <- data.frame(LValue = character(),
                              SampleSize = numeric(),
                              Mean = numeric(),
                              SD = numeric(),
                              Variance = numeric(),
                              Skewness = numeric(),
                              Kurtosis = numeric(),
                              adpvalue = numeric())
  
  for (s in sample.size) {
    TestStat1 <- numeric(R)
    
    # Simulaciones para datos homogéneos
    for (r in 1:R) {
      z <- gamma_sar_sample(L, mu, s)
      TestStat1[r] <- ( sd(z)/mean(z) ) - (1 / sqrt(L))#S_CV(z, L) # Calcular S_CV
    }
    
    # Estimar E[S_CV] y Var[S_CV]
    mean_S_CV <- mean(TestStat1)
    var_S_CV <- var(TestStat1)
    
    # Calcular Z_CV normalizado
    Z_CV <- (TestStat1 - mean_S_CV) / sqrt(var_S_CV)
    
    TestStatistics1[[as.character(s)]] <- data.frame("SampleSize" = rep(s, R), "Test_Statistics" = Z_CV)
    
    # Estadísticas resumen
    mean_val <- mean(Z_CV)
    sd_val <- sd(Z_CV)
    var_val <- var(Z_CV)
    skewness_val <- skewness(Z_CV)
    kurtosis_val <- kurtosis(Z_CV)
    ad_p_value <- ad.test(TestStatistics1[[as.character(s)]]$Test_Statistics)$p.value
    
    summary_stats <- rbind(summary_stats, data.frame(LValue = as.character(L),
                                                     SampleSize = s,
                                                     Mean = mean_val,
                                                     SD = sd_val,
                                                     Variance = var_val,
                                                     Skewness = skewness_val,
                                                     Kurtosis = kurtosis_val,
                                                     adpvalue = ad_p_value))
  }
  
  all_TestStatistics[[as.character(L)]] <- TestStatistics1
  all_summary_stats[[as.character(L)]] <- summary_stats
  
  save(all_TestStatistics, all_summary_stats, file = paste0("./Data/resultsZCV2_", L, ".Rdata"))
}
