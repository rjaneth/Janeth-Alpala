library(png)
source("../imagematrix.R")



#reparametrized distribution for mu
rGI0 <- function(n, p_alpha, p_mu, p_Looks) {
  return(rgamma(n, shape = p_Looks, rate = p_Looks) / rgamma(n, shape = -p_alpha, rate = p_mu * (-p_alpha - 1)))
}


Z <- Phantom1 <- readPNG("../../../Figures/PNG/Phantom1.png")[,,1]
dim(Phantom1)
#plot(imagematrix(Phantom1))

E <- function(a, g) {-g/(a+1)}
H <- function(a, g) {(a-1)/a - log(-a/g, base=exp(1))}
HE1 <- function(a) {(a-1)/a - log(a/(a+1), base=exp(1))}

HE1(-10)
H(-2,50)
H(-2,5)
E(-6, 10)

# Phantom

## Return
## z: return
## up: upper, lo: lower
## le: left, ri: right
## b: black, w: white

p.up.le <- Phantom1[1:250,1:250]
p.up.ri <- Phantom1[1:250,251:500]
p.bo.le <- Phantom1[251:500,1:250]
p.bo.ri <- Phantom1[251:500,251:500]


### Upper Left (means 1 and 1)
z.up.le.b <- rGI0(n=sum(1-p.up.le), 
                  p_alpha=-2, p_mu=1, p_Looks=5)
mean(z.up.le.b)

z.up.le.w <- rGI0(n=sum(p.up.le), 
                  p_alpha=-10, p_mu=9, p_Looks=5)
mean(z.up.le.w)
z.up.le <- p.up.le
z.up.le[p.up.le==0] <- z.up.le.b
z.up.le[p.up.le==1] <- z.up.le.w
plot(imagematrix(equalize(z.up.le)))

### Upper Right (means 50 and 5)
z.up.ri.b <- rGI0(n=sum(1-p.up.ri), 
                  p_alpha=-10, p_mu=1, p_Looks=5)
mean(z.up.ri.b)

z.up.ri.w <- rGI0(n=sum(p.up.ri), 
                  p_alpha=-2, p_mu=9, p_Looks=5)
mean(z.up.ri.w)
z.up.ri <- p.up.ri
z.up.ri[p.up.ri==0] <- z.up.ri.b
z.up.ri[p.up.ri==1] <- z.up.ri.w
plot(imagematrix(equalize(z.up.ri)))


### Bottom Right (means 10 and 10)

z.bo.ri.b <- rGI0(n=sum(1-p.bo.ri), 
                  p_alpha=-9, p_mu=1, p_Looks=5)
mean(z.bo.ri.b)

z.bo.ri.w <- rGI0(n=sum(p.bo.ri), 
                  p_alpha=-3, p_mu=9, p_Looks=5)
mean(z.bo.ri.w)
z.bo.ri <- p.bo.ri
z.bo.ri[p.bo.ri==0] <- z.bo.ri.b
z.bo.ri[p.bo.ri==1] <- z.bo.ri.w
plot(imagematrix(equalize(z.bo.ri)))

### Bottom left (means 10 and 10)

z.bo.le.b <- rGI0(n=sum(1-p.bo.le), 
                  p_alpha=-3, p_mu=1, p_Looks=5)
mean(z.bo.le.b)
z.bo.le.w <- rGI0(n=sum(p.bo.le), 
                  p_alpha=-9, p_mu=9, p_Looks=5)
mean(z.bo.le.w)
z.bo.le <- p.bo.le
z.bo.le[p.bo.le==0] <- z.bo.le.b
z.bo.le[p.bo.le==1] <- z.bo.le.w
plot(imagematrix(equalize(z.bo.le)))

Z[1:250,1:250] <- z.up.le
Z[1:250,251:500] <- z.up.ri
Z[251:500,1:250] <- z.bo.le
Z[251:500,251:500] <- z.bo.ri

plot(imagematrix(equalize(Z)))



save(Z, file="./Data/Phantom_4_regions.Rdata")


