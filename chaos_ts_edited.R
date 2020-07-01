
library(ineq)

data(Ilocos)
attach(Ilocos)
## extract and rescale income for the provinces "Pangasinan" und "La Union"
income.p <- income[province=="Pangasinan"]/10000
income.u <- income[province=="La Union"]/10000
## compute the Lorenz curves
Lc.p <- Lc(income.p)
Lc.u <- Lc(income.u)

## plot both Lorenz curves

plot(Lc.p)
lines(Lc.u, col=2)

data("AirPassengers")

plot(Lc(AirPassengers))


## above scrape ##

# this works

library(fNonlinear)

lorentzSim(times = seq(0, 20, by = 0.01), parms = c(sigma = 16, r = 45.92, b = 4), 
           start = c(-14, -13, 47), doplot = TRUE)

lorentzSim(times = seq(0, 40, by = 0.01), parms = c(sigma = 16, r = 45.92, b = 4), 
           start = c(1, 2, 3), doplot = TRUE)

lorenz <- lorentzSim(
  times = seq(0, 40, by = 0.01),
  parms = c(sigma = 16, r = 45.92, b = 4),
  start = c(-14, -13, 47),
  doplot = TRUE)

recurrencePlot(lorentz[, 2], m = 3, d = 2, end.time = 800, eps = 3,
               nt = 5, pch = '.', cex = 2, radius = 0.5)

head(lorenz) # t, x, y, z useful for creating verticle or horizontal plots

plot(lorenz[, 2], lorenz[, 3]) # base methods


## to create ts data

x <- read.csv(file.choose())
names(x)

x <- abs(round(rnorm(100), 2))
myts <- ts(x) #, start=c(2017, 07), end=c(2017, 24)) 
myts
plot(myts)

# explorative analysis

x <- log10(myts)

#par(mfrow = c(2, 1)) #mar = c(0, 0, 0, 0))
#par()

plot(x, ax = F)
box()
plot(x[length(x):1], type = "l", ax = F)
box()


par(mfrow = c(2, 1), mar = c(2, 2, 0, 0))

autopairs(x, lag = 1, type = "regression")
autopairs(x, lag = 3, type = "regression")


library(tseriesChaos)

output <-lyap_k(x, m=3, d=2, s=00, t=40, ref=1700, k=2, eps=4)
plot(output)
lyap(output, 0.73, 2.47)


## ploting lorenz attractor nonlineartseries

suppressMessages(library('nonlinearTseries'))

library('plot3D')

lorenz(sigma = 10, beta = 8/3, rho = 28, start = c(-13, -14, 47),
       time = seq(0, 50, by = 0.01), do.plot = TRUE)

# let's plot the phase space of the simulated lorenz system
scatter3D(lor$x, lor$y, lor$z,
          main = "Lorenz's system phase space",
          col = 1, type="o",cex = 0.3)


maxLyapunov(myts, min.embedding.dim = 2,
            max.embedding.dim = 2, time.lag = 1, radius,
            theiler.window = 1, min.neighs = 5, min.ref.points = 500,
            max.time.steps = 10, number.boxes = NULL, sampling.period = 1,
            do.plot = TRUE)


lorenz = function(t, x, parms = c(sigma = 10, b = 8/3, r = 28)) {
  
  X = x[1] # not correct
  Y = x[2] # not correct
  Z = x[3] # not correct
  
  with(as.list(parms), {
    dX = parms["sigma"] * (Y - X)
    dY = -x * Z + parms["r"] * X - Y
    dZ = X * Y - parms["b"] * Z
    out <- list(c(dX, dY, dZ))
    
  })
  plot(out[2], out[3])
}

lorenz(1:100, ts(rnorm(100)))


lyap_exp <- function(x, l){
  y <- diff(x, l)
  y <- abs(y)
  
  lamb.e <- matrix(NA, length(y))
  for(i in 1:length(y)){
    lamb.e[i] <- (1/length(y)*log(y[i+1]/y[i]))  
  }
  
  lamb.v = as.vector(lamb.e)
  mle <- max(lamb.e)
  
  # kolmogorov-sinai entropy
  
  kse <- matrix(NA, length(lamb.v))
  for(i in 1:length(lamb.v)){
    if(lamb.v[i] > 0) kse[i]=lamb.v[i] 
  }
  
  #kse.sum <- sum(kse, na.rm = TRUE)
  
  scatter.smooth(lamb.e)
  
  ans = list(exponents=lamb.v, mle=mle)
  return(ans)
}

# edited version of lyapunov's exponent

lyap_exp <- function(x, l){
  
  y <- diff(x, l)
  y <- abs(y)
  y <- ts(y)
  
  lamb.e <- matrix(NA, length(y))
  
  for(i in 1:length(y)){
    lamb.e[i] <- (1/length(y)*log(y[i+1]/y[i]))  
  }
  
  lamb.v = as.vector(na.omit(lamb.e))
  mle <- max(lamb.v)
  
  # kolmogorov-sinai entropy
  
  kse <- matrix(length(lamb.v))
  for(i in 1:length(lamb.v)){
    if(lamb.e[i] > 0){kse[i]=lamb.v[i]}  
  }
  
  kse.sum <- sum(kse, na.rm = TRUE)
  
  # plot the lyapunov es'
  
  scatter.smooth(lamb.e)
  
  # return
  
  ans = list(exponents=lamb.v, mle=mle, divergence=kse.sum)
  return(ans)
  
}