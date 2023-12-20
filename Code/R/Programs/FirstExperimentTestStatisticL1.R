library(ggplot2)
library(ggsci)

source("../MainFunctions/ebrahimi_estimator.r")
source("../MainFunctions/gamma_sar_sample.r")




set.seed(1234567890, kind="Mersenne-Twister")

R <- 3000 
mu <- 1
sample.size <- c(9, 25, 49, 121)

TestStatistics <- NULL

for(s in sample.size){
  for(r in 1:R){
    
    z <- gamma_sar_sample(1, mu, s)
    TestStat <- ebrahimi_estimator(z) - 1 - log(mean(z))
    TestStatistics <- rbind(TestStatistics,
                            c(s, TestStat))
  }
}

TestStatistics <- data.frame(TestStatistics)
names(TestStatistics) <- c("Sample Size", "Test Statistics")

# ggplot(TestStatistics, aes(x=`Test Statistics`, col=`Sample Size`, group=`Sample Size`)) +
#   geom_density()


ggplot(TestStatistics, aes(x = `Test Statistics`, col = factor(`Sample Size`), linetype = factor(`Sample Size`))) +
  geom_line(stat = "density", linewidth = 1.0) +  
  scale_color_manual(
    values = pal_cosmic()(length(sample.size)),
    name = "Sample Size"
  ) +
  scale_linetype_manual(
    values = rep("solid", length(sample.size)),  
    name = "Sample Size"
  ) +
  labs(x = "Test Statistics", y = "Density") +
  theme_minimal() +
  theme(
    text = element_text(family = "serif"),
    legend.position = "bottom",
    # legend.title = element_text(size = 10),
    # legend.text = element_text(size = 8),
    legend.key.size = unit(1, "lines")  
  )
