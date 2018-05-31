# DynGamma.R - functions to compute gamma_u, gamma_v, and their
#              first derivates for the dynamic model
#

Gamma_nu <- function(rho,T) {
  Gamma_nu <- d.Gamma_nu <- matrix(0,nrow=T, ncol=T) 
  if (rho==0) return(list(Gamma_nu= Gamma_nu, d.Gamma_nu = d.Gamma_nu)) 
  diag(Gamma_nu)      <- cumsum(rho^(2*(0:(T-1))))
  diag(d.Gamma_nu)    <- cumsum((2*(0:(T-1)))*rho^((2*(0:(T-1))-1)))
  for (j in 2:T) {
    for (i in 1:j-1) {
      Gamma_nu[i,j]   <- (rho^(j-i))*sum(rho^(2*(0:(i-1))))
      Gamma_nu[j,i]   <- Gamma_nu[i,j]
      d.Gamma_nu[i,j] <- (rho^(j-i-1))*sum((j-i+2*(0:(i-1)))*rho^(2*(0:(i-1))))
      d.Gamma_nu[j,i] <- d.Gamma_nu[i,j]
    }
  }
  return(list(Gamma_nu= Gamma_nu, d.Gamma_nu = d.Gamma_nu))
}
  
Gamma_u <- function(rho, T, model="rao") {
  if (model=="rao") {
    Gamma_u   <- diag(T)
    d.Gamma_u <- 0*diag(T)
    return(list(Gamma_u=Gamma_u,d.Gamma_u=d.Gamma_u))
  } else if (model=="dyn") {
    rho_power <- rho^(0:(T-1))
    Gamma_u   <- rho_power %*% t(rho_power)
    rho_d     <- c(0,(rho^(0:(T-2)))*(1:(T-1)))
    d.Gamma_u <- rho_power %*% t(rho_d) + rho_d %*% t(rho_power)
    return(list(Gamma_u=Gamma_u,d.Gamma_u=d.Gamma_u))
  } else {
    print("Must indicate appropriate model rao or dyn") 
  }
}
