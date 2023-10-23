
rm(list = ls())
library(gtools)
library(stats4)
library(rmutil)
library(gsl)
library(invgamma)
library(ggplot2)
library(ggthemes)
theme_set(theme_pander() +
            theme(text=element_text(family="serif"),
                  legend.position = "top")
)


alpha <- -2
mu <- 1
L <- 1

gam=-mu*(alpha+1)

r <- 500
n_values <- c(9, 25, 49, 81, 200,1000, 5000,  10000, 20000) 




means_entropies <- c()
means_entropies2 <- c()
entropies2 <- c()

for (n in n_values) {
  entropies <- c()
  entropies2 <- c()
  for (i in 1:r) {
    sample_data <- gamma_sar_sample(L, mu, n)
    GI0_samples <- Z_samples(L, alpha, mu,  n)
    Np_entropy_gamma <- van_es_estimator(sample_data)
    NP_entropy_GI0 <- van_es_estimator(GI0_samples)
    entropies <- c(entropies, Np_entropy_gamma)
    entropies2 <- c(entropies2, NP_entropy_GI0)
    log_gamma_sar<- log(dgamma(sample_data, shape = L, rate = L / mu))
    #log_var <- var(log_sample_data)
    
  }
  #Gamma SAR
  mean_entropy <- mean(entropies)
  means_entropies <- c(means_entropies, mean_entropy)
  # GIO
  mean_entropy2 <- mean(entropies2)
  means_entropies2 <- c(means_entropies2, mean_entropy2)
}


# plot(n_values, means_entropies, type = "o", main = "Convergence of Mean Non-Parametric Entropies",
     # xlab = "Sample Size (n)", ylab = "Mean Entropy")
# grid()

### Using ggplot2

df.MeanGamma <- data.frame(n_values, means_entropies)
ggplot(df.MeanGamma, aes(x=n_values, y=means_entropies)) +
  geom_hline(yintercept = entropy_gamma_sar(1, 1)) +
  geom_line() +
  geom_point() +
  xlab("Sample Size") +
  scale_x_continuous(breaks=n_values) +
  ylab("Mean Non-parametric Entropy")
  

entropy_Gamma_sar <- entropy_gamma_sar(L, mu)
cat("mean of Non-Parametric Entropies gamma sar :", means_entropies, "\n")
cat("Analitical Entropy gamma sar:", entropy_Gamma_sar, "\n")

entropy_GI0 <- entropy_gI0(mu, alpha, L)
cat("Analitical Entropy GIO:", entropy_GI0, "\n")

cat("mean of Non-Parametric Entropies GI0 :", means_entropies2, "\n")

second <- second_term(alpha, L)
cat("second term of GI0 :", second, "\n")



# # test
#   alpha1 <- 0.05
#   D0 <- entropy_Gamma_sar
# 
# 
#   Z <- (means_entropies - D0) / (sqrt(var(log_gamma_sar)) / sqrt(n_values))
# 
#   # Valor crítico
#   z_alpha_2 <- qnorm(1 - alpha1/2)
# 
# 
# 
# 
#   p_values <- 2 * (1 - pnorm(abs(Z)))
# 
# 
#   rechazo <- p_values < alpha1
# 
#   # Resultados
#   resultados <- data.frame(
#     Sample_Size = n_values,
#     Mean_Entropy = means_entropies,
#     Z_Statistic = Z,
#     P_Value = p_values,
#     Rechazo_H0 = rechazo
#   )
# 
#   print(resultados)
# 
