difference.betweeen.GammaSAR.GI0.n.Look <- function(alpha, L) {
  return(
    -L - lgamma(L-alpha) + 
      (L-alpha)*(digamma(L - alpha))- 
      (1-alpha)*digamma(- alpha)+
      log(-1 - alpha)+
      lgamma(-alpha)
  )
}

x=difference.betweeen.GammaSAR.GI0.n.Look(-1.5, 2)