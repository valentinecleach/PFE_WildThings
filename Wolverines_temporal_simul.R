
## Initialisation

omega <- array(0, dim = c(4, 4, n.years1))

## Boucle temporelle
for (t in 1:n.years1) {
  
  gamma[t] <- runif(1, 0, 1)
  w[t]     <- runif(1, 0, 0.59)
  h[t]     <- runif(1, 0, 0.39)
  phi[t]   <- 1 - h[t] - w[t]
  
  omega[1, , t] <- c(1 - gamma[t], gamma[t], 0, 0)
  omega[2, , t] <- c(0, phi[t], h[t], w[t])
  omega[3, , t] <- c(0, 0, 0, 1)
  omega[4, , t] <- c(0, 0, 0, 1)
}

omega
