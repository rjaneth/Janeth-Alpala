library(plotly)
library(latex2exp)

# Define la función de entropía
entropy_gI0 <- function(mu, alpha, L) {
  term1 <- L - log(L) + log(gamma(L)) + (1 - L) * digamma(L) + log(mu)   
  term2 <- -L - log(gamma(L-alpha)) + (L-alpha)*(digamma(L - alpha))- (1-alpha)*digamma(- alpha)+log(-1 - alpha)+log(gamma(-alpha))
  entropy <- term1 + term2 
  return(entropy)
}

mu <- seq(1, 20, length.out = 200)
alpha <- seq(-20, -2, length.out = 200)
# entropyL1 <- outer(mu, alpha, function(mu, alpha) entropy_gI0(mu, alpha, L=1))
# entropyL3 <- outer(mu, alpha,  function(mu, alpha) entropy_gI0(mu, alpha, L=3))
# entropyL8 <- outer(mu, alpha,  function(mu, alpha) entropy_gI0(mu, alpha, L=8))
# entropyL12 <- outer(mu, alpha, function(mu, alpha) entropy_gI0(mu, alpha, L=12))
# entropyL100 <- outer(mu, alpha, function(mu, alpha) entropy_gI0(mu, alpha, L=100))

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



fig <- fig %>% layout(
  scene = list(
    xaxis = list(title = "\u03BC"), # \u03BC representa el símbolo μ en Unicode
    yaxis = list(title = "\u03B1"), # \u03B1 representa el símbolo α en Unicode
    zaxis = list(title = "Entropy"),
    annotations = list(
      list(
        text = "L=100", x = 11, y = -23,
        showarrow = FALSE,
        xshift = -15,
        xanchor = "center",
        bgcolor = ""
      ),
      list(
        text = "L=12", x = 11, y = -17,
        showarrow = FALSE,
        xshift = -15,
        xanchor = "center",
        bgcolor = ""
      ),
      list(
        text = "L=8", x = 12, y = -10,
        showarrow = FALSE,
        xshift = -15,
        xanchor = "center",
        bgcolor = ""
      ),
      list(
        text = "L=3", x = 13, y = -5,
        showarrow = FALSE,
        xshift = -15,
        xanchor = "center",
        bgcolor = ""
      ),
      list(
        text = "L=1", x = 14, y = 5,
        showarrow = FALSE,
        xshift = -15,
        xanchor = "center",
        bgcolor = ""
      )
    )
  )
)

camera_position <- list(
  eye = list(x = -2, y = 0, z = 0),
  center = list(x = -2, y = 0, z = 0),
  up = list(x = 0, y = 0, z = 1 )
)
# Muestra el gráfico
fig
orca(fig, file = "mi_grafico3.pdf")
