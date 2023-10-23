# This code runs examples of the sample mean of the non-parametric...

set.seed(1234567890, kind="Mersenne-Twister")

sample_sizes <- c(9, 25, 121, 1000, 10000, 100000, 1000000)
R <- 500
mu <- 1
L <- 1

output <- NULL

for(ssize in sample_sizes) {
  v.nonparametric.entropy <- NULL
  for(r in 1:R) {
    sample <- gamma_sar_sample(L, mu, ssize)
    v.nonparametric.entropy[r] <- van_es_estimator(sample)
  }
  output <- rbind(output, c(mean(v.nonparametric.entropy), sd(v.nonparametric.entropy)))
}

output <- data.frame(output)
names(output) <- c("Mean", "Std")

output$SampleSize <- sample_sizes

ggplot(output, aes(x=SampleSize, y=Mean)) +
  geom_hline(yintercept = entropy_gamma_sar(1, 1)) +
  geom_line() +
  geom_point() +
  xlab("Sample Size") +
  scale_x_log10(breaks=output$SampleSize) +
  ylab("Mean Non-parametric Entropy")


### Show curves for same mu=1 and L in 1, 3, 8, 12.
### Add confidence intervals to each point based on the standard deviation
### Save functions on a different file; use "standard R names"
### Each function that produces a plot must be in a different file
