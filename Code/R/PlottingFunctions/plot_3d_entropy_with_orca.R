library(plotly)
library(latex2exp)

# Entropy
entropy_gI0 <- function(mu, alpha, L) {
  term1 <- L - log(L) + log(gamma(L)) + (1 - L) * digamma(L) + log(mu)   
  term2 <- -L - log(gamma(L-alpha)) + (L-alpha)*(digamma(L - alpha))- (1-alpha)*digamma(- alpha)+log(-1 - alpha)+log(gamma(-alpha))
  entropy <- term1 + term2 
  return(entropy)
}

mu <- seq(1, 10, length.out = 200)
alpha <- seq(-20, -2, length.out = 200)

entropyL1 <- outer(alpha, mu, function(alpha, mu) entropy_gI0(mu, alpha, L = 1))
entropyL3 <- outer(alpha, mu, function(alpha, mu) entropy_gI0(mu, alpha, L = 3))
entropyL8 <- outer(alpha, mu, function(alpha, mu) entropy_gI0(mu, alpha, L = 8))
entropyL12 <- outer(alpha, mu, function(alpha, mu) entropy_gI0(mu, alpha, L = 12))
entropyL100 <- outer(alpha, mu, function(alpha, mu) entropy_gI0(mu, alpha, L = 100))

fig <- plot_ly(x = mu, y = alpha, showscale = FALSE)

fig <- fig %>% add_surface(z = ~entropyL1, colorscale = list(c(0, 1), c("rgb(255,107,184)", "rgb(128,0,64)")), name = "L=1")
fig <- fig %>% add_surface(z = ~entropyL3, colorscale = list(c(0, 1), c("rgb(247,226,157)", "rgb(255,192,7)")), name = "L=3")
fig <- fig %>% add_surface(z = ~entropyL8,colorscale = list(c(0, 1), c("rgb(107,255,184)", "rgb(0,124,90)")),  name = "L=8")
fig <- fig %>% add_surface(z = ~entropyL12,colorscale = list(c(0,1),c("rgb(129,212,247)","rgb(12,177,247)")),  name = "L=12")
fig <- fig %>% add_surface(z = ~entropyL100, colorscale = list(c(0, 1), c("rgb(182,142,242)", "rgb(104,3,255)")), name = "L=100")

# 
fig <- fig %>% layout(
  scene = list(
    camera = list(
      eye = list(x = -1.5, y = -1.5, z = 1.5),  
      center = list(x = 0, y = 0, z = 0),  
      up = list(x = 0, y = 0, z = 1)  
    ),
    xaxis = list(title = "\u03BC"),
    yaxis = list(title = "\u03B1"),
    zaxis = list(title = "Entropy")
  )
)
fig

#  PDF
orca(fig, file = "Entropy_plot_3d_1.pdf")