rm(list = ls())

library(ggplot2)

library(grid)
library(ggsci)
library(invgamma)
library(latex2exp)
library(gridExtra)
source("../MainFunctions/ebrahimi_estimator.R")
source("../MainFunctions/gamma_sar_sample.R")
source("../MainFunctions/correa_estimator.R")
source("../MainFunctions/bootstrap_correa_estimator.R")
source("../MainFunctions/bootstrap_correa_estimator_log_mean.R")

source("../MainFunctions/al_omari_1_estimator.R")
source("../MainFunctions/bootstrap_al_omari_1_estimator.R")

source("../MainFunctions/entropy_gI0.R")
source("../MainFunctions/gi0_sample.R")

set.seed(1234567890, kind = "Mersenne-Twister")

R <- 3000
mu <- 1
L <- 5
B <- 2
sample.size <- c( 81)

alpha1_values <- c( -3, -8, -20, -100)


plots_list <- list()

for (alpha_val in alpha1_values) {

  alpha1 <- alpha_val
  
  TestStatistics1 <- NULL
  TestStatistics2 <- NULL
  
  mean_entropy1 <- numeric(length(sample.size))
  mean_entropy2 <- numeric(length(sample.size))
  
  for (s in sample.size) {
    TestStat1 <- numeric(R)
    TestStat2 <- numeric(R)
    
    for (r in 1:R) {
      z <- gi0_sample(mu, alpha1, L, s)
      
      TestStat1[r] <- bootstrap_correa_estimator_log_mean(z, B) - (L - log(L) + lgamma(L) + (1 - L) * digamma(L))
      TestStat2[r] <- bootstrap_correa_estimator_log_mean(z, B) - (L - log(L) + lgamma(L) + (1 - L) * digamma(L) - L - lgamma(L - alpha1) + (L - alpha1) * (digamma(L - alpha1))- (1 - alpha1) * digamma(-alpha1) + log(-1 - alpha1) + lgamma(-alpha1))
      
      #TestStat2[r] <- -bootstrap_correa_estimator_log_mean(z, B) - (-L + log(L) - lgamma(L) - (1 - L) * digamma(L) + L + lgamma(L - alpha1) - (L - alpha1) * (digamma(L - alpha1)) + (1 - alpha1) * digamma(-alpha1) - log(-1 - alpha1) - lgamma(-alpha1))
    }
    
    mean_entropy1[sample.size == s] <- mean(TestStat1)
    cat(" mean for TestStat1", mean_entropy1[sample.size == s], "\n")
    
    mean_entropy2[sample.size == s] <- mean(TestStat2)
    cat(" mean for TestStat2", mean_entropy2[sample.size == s], "\n")
    
    TestStatistics1 <- rbind(TestStatistics1, data.frame("Sample_Size" = rep(s, R), "Test_Statistics" = TestStat1))
    TestStatistics2 <- rbind(TestStatistics2, data.frame("Sample_Size2" = rep(s, R), "Test_Statistics2" = TestStat2))
  }
  
  difference.betweeen.GammaSAR.GI0.n.Look <- function(alpha, L) {
    return(
      -L - lgamma(L-alpha) + 
        (L-alpha)*(digamma(L - alpha))- 
        (1-alpha)*digamma(- alpha)+
        log(-1 - alpha)+
        lgamma(-alpha)
    )
  }
  
  max_density1 <- max(density(TestStatistics1$Test_Statistics)$y, na.rm = TRUE)
  max_density2 <- max(density(TestStatistics2$Test_Statistics2)$y, na.rm = TRUE)
  
  plot <- ggplot(TestStatistics1, aes(x = Test_Statistics)) +
    geom_line(aes(col = factor(Sample_Size), linetype = factor(Sample_Size)), stat = "density", linewidth = 1.0) +
    geom_line(data = TestStatistics2, aes(x = Test_Statistics2), stat = "density", color = "blue", linetype = "dashed", linewidth = 1) +
    # annotate("text", x = max(TestStatistics1$Test_Statistics) + 0.1, y = max_density1,
    #          label = TeX("$\\Gamma_{\\tiny{SAR}}$"), vjust = 5, hjust = 5, color = "red") +
    # annotate("text", x = max(TestStatistics2$Test_Statistics2) + 0.1, y = max_density2,
    #          label = TeX("$G_I^0$"), vjust = 5, hjust = 5, color = "blue") +
    geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.n.Look(-3, 5), col="darkgreen", linetype = "dotted", linewidth = 0.6) +
    annotate("text", x=difference.betweeen.GammaSAR.GI0.n.Look(-3, 5), y=7, parse=TRUE, label="alpha==-3", hjust=-.1, vjust=-0.5, col="darkgreen") +
    geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.n.Look(-8, 5), col="darkgreen", linetype = "dotted", linewidth = 0.6) +
    annotate("text", x=difference.betweeen.GammaSAR.GI0.n.Look(-8, 5), y=7, parse=TRUE, label="alpha==-8", hjust=-.1, vjust=1.0, col="darkgreen") +
    geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.n.Look(-20, 5), col="darkgreen", linetype = "dotted", linewidth = 0.6) +
    annotate("text", x=difference.betweeen.GammaSAR.GI0.n.Look(-20, 5), y=7, parse=TRUE, label="alpha==-20", hjust=-.1, vjust=2.0, col="darkgreen") +
    geom_vline(xintercept = difference.betweeen.GammaSAR.GI0.n.Look(-1000, 5), col="darkgreen", linetype = "dotted", linewidth = 0.6) +
    annotate("text", x=difference.betweeen.GammaSAR.GI0.n.Look(-1000, 5), y=7, parse=TRUE, label="alpha==-1000", hjust=-.1, vjust=4.0, col="darkgreen") +
    labs(x = "Test Statistics", y = "Empirical Density") +
    #annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("alpha == %s", alpha_val)), hjust = 1.2, vjust = 1.02, size = 3) +
    theme_minimal() +
    theme(
      text = element_text(family = "serif"),
      legend.position = "none"  
    )
  
  

  plots_list[[as.character(alpha_val)]] <- plot
}

grid.arrange(grobs = plots_list, ncol = 2)

# plot <- ggplot(TestStatistics1, aes(x = Test_Statistics, col = factor(Sample_Size), linetype = factor(Sample_Size))) +
#   geom_line(stat = "density", linewidth = 1.0) +
#   scale_color_manual(values = rainbow(7)[1:2], name = "Sample Size") +
#   scale_linetype_manual(values = rep("solid", length(sample.size)), name = "Sample Size") +
#   geom_line(data = TestStatistics2, aes(x = Test_Statistics2, col = factor(Sample_Size2), linetype = factor(Sample_Size2)), stat = "density", color = "blue", linetype = "dashed", linewidth = 1) +
#   annotate("text", x = max(TestStatistics1$Test_Statistics) + 0.1, y = max_density1,
#            label = TeX("$\\Gamma_{\\tiny{SAR}}$"), vjust = 5, hjust = 5, color = "red") +
#   annotate("text", x = max(TestStatistics2$Test_Statistics2) + 0.1, y = max_density2,
#            label = TeX("$G_I^0$"), vjust = 5, hjust = 5, color = "blue") +
#   labs(x = "Test Statistics", y = "Density") +
#   annotate("text", x = Inf, y = Inf, label = parse(text = sprintf("alpha == %s", alpha_val)), hjust = 1.2, vjust = 1.02, size = 3) +
#   theme_minimal() +
#   theme(
#     text = element_text(family = "serif"),
#     legend.position = "bottom",
#     legend.key.size = unit(1, "lines")
#   ) 
  
#   # Crear gráfico y agregar a la lista
#   plot <- ggplot(TestStatistics1, aes(x = Test_Statistics, col = factor(Sample_Size), linetype = factor(Sample_Size))) +
#     geom_line(stat = "density", linewidth = 1.0) +
#     scale_color_manual(values = rainbow(7)[1:2], name = "Sample Size") +
#     scale_linetype_manual(values = rep("solid", length(sample.size)), name = "Sample Size") +
#     geom_line(data = TestStatistics2, aes(x = Test_Statistics2, col = factor(Sample_Size2), linetype = factor(Sample_Size2)), stat = "density", color = "blue", linetype = "dashed", linewidth = 1) +
#     annotate("text", x = max(TestStatistics1$Test_Statistics) + 0.1, y = max_density1,
#              label = TeX("$\\Gamma_{\\tiny{SAR}}$"), vjust = 5, hjust = 5, color = "red") +
#     annotate("text", x = max(TestStatistics2$Test_Statistics2) + 0.1, y = max_density2,
#              label = TeX("$G_I^0$"), vjust = 5, hjust = 5, color = "blue") +
#     labs(x = "Test Statistics", y = "Density") +
#     theme_minimal() +
#     theme(
#       text = element_text(family = "serif"),
#       legend.position = "bottom",
#       legend.key.size = unit(1, "lines")
#     )
#   
#   # Agregar el gráfico a la lista
#   plots_list[[as.character(alpha_val)]] <- plot
# }
# 
# # Mostrar los subgráficos
# 
# grid.arrange(grobs = plots_list, ncol = 2)

#   # Gráficos para cada valor de alpha1
#   max_density1 <- max(density(TestStatistics1$Test_Statistics)$y, na.rm = TRUE)
#   max_density2 <- max(density(TestStatistics2$Test_Statistics2)$y, na.rm = TRUE)
#   
#   ggplot(TestStatistics1, aes(x = Test_Statistics, col = factor(Sample_Size), linetype = factor(Sample_Size))) +
#     geom_line(stat = "density", linewidth = 1.0) +
#     scale_color_manual(values = rainbow(7)[1:2], name = "Sample Size") +
#     scale_linetype_manual(values = rep("solid", length(sample.size)), name = "Sample Size") +
#     geom_line(data = TestStatistics2, aes(x = Test_Statistics2, col = factor(Sample_Size2), linetype = factor(Sample_Size2)), stat = "density", color = "blue", linetype = "dashed", linewidth = 1) +
#     geom_line(data = TestStatistics3, aes(x = Test_Statistics3, col = factor(Sample_Size3), linetype = factor(Sample_Size3)), stat = "density", color = "green", linetype = "dashed", linewidth = 1) +
#     geom_line(data = TestStatistics4, aes(x = Test_Statistics4, col = factor(Sample_Size4), linetype = factor(Sample_Size4)), stat = "density", color = "black", linetype = "dashed", linewidth = 1) +
#     annotate("text", x = max(TestStatistics1$Test_Statistics) + 0.1, y = max_density1,
#              label = TeX("$\\Gamma_{\\tiny{SAR}}$"), vjust = 5, hjust = 5, color = "red") +
#     annotate("text", x = max(TestStatistics2$Test_Statistics2) + 0.1, y = max_density2,
#              label = TeX("$G_I^0$"), vjust = 5, hjust = 5, color = "blue") +
#     labs(x = "Test Statistics", y = "Density") +
#     theme_minimal() +
#     theme(
#       text = element_text(family = "serif"),
#       legend.position = "bottom",
#       legend.key.size = unit(1, "lines")
#     )
#   
#   # Guardar o mostrar el gráfico aquí si es necesario
# }
