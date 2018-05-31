######################################################################################################### 

# THIS PROGRAM CONTAINS THE FUNCTIONS USED IN THE Rao-Yu/Dynamic 

# Autheur: Mamadou S. Diallo
# Date: January 2014

######################################################################################################## 


### This will load the library where are stored the packages .libPaths('C:/Program
### Files/R/Cont_Library')

### Loading the libraries needed for this program ###

library(Matrix)
library(MASS)



### Remove all objects to start fresh rm(list=ls())


## Theta_hat - generate the direct estimates for the dynamic and Rao-Yu model
covar <- function(rho_cov, var1, var2, M, T) {
    vcov_e <- 0 * diag(n)
    Gamma <- diag(T)
    for (t in 1:(T - 1)) {
        Gamma[t, (t + 1):T] <- rho_cov^(1:(T - t))
        Gamma[(t + 1):T, t] <- Gamma[t, (t + 1):T]
    }
    var <- seq(var1, var2, length.out = n)
    for (d in 1:M) {
        sqrt_var_vec <- as.vector(sqrt(var[(1 + (d - 1) * T):(d * T)]))
        vcov_e[(1 + (d - 1) * T):(d * T), (1 + (d - 1) * T):(d * T)] <- Gamma * ((sqrt_var_vec) %*% 
            t(sqrt_var_vec))
    }
    diag(vcov_e) <- var
    return(vcov_e)
}

## X <- cbind(rep(1,n), rbinom(n,1,p1=0.2), rbinom(n,1,p2=0.8)) #generate the auxiliary variables

Theta_hat <- function(sig2_nu, sig2_u, rho, T, M, beeta, X, vcov_e, Stationarity) {
    n <- T * M
    u_rao <- rep(rnorm(M, mean = 0, sd = sqrt(sig2_u)), each = T)
    u_dyn <- rho^(T - 1) * rep(rnorm(M, mean = 0, sd = sqrt(sig2_u)), each = T)
    nu <- matrix(0, M, T)
    if (toupper(Stationarity) == "FALSE") {
        nu[, 1] <- rnorm(M, mean = 0, sd = sqrt(sig2_nu))
        for (t in 2:T) nu[, t] <- rho * nu[, t - 1] + rnorm(M, mean = 0, sd = sqrt(sig2_nu))
    } else if (toupper(Stationarity) == "TRUE") {
        # arima.sim produces a series of variance 1 so it is necessary to multiply by sqrt(sig2_nu) to
        # get the right variance
        for (m in 1:M) nu[m, ] <- arima.sim(n = T, model = list(order = c(1, 0, 0), ar = rho, sd = sqrt(sig2_nu))) * 
            sqrt(sig2_nu)
    } else {
        print("Error: must indicate stationarity status (TRUE or FALSE)")
        break
    }
    nu <- as.vector(t(nu))
    theta_rao <- as.vector(X %*% beeta + u_rao + nu)
    theta_dyn <- as.vector(X %*% beeta + u_dyn + nu)
    e <- mvrnorm(n = 1, mu = rep(0, n), Sigma = vcov_e)
    y_hat_rao <- as.vector(theta_rao + e)
    y_hat_dyn <- as.vector(theta_dyn + e)
    cv_rao <- sqrt(diag(vcov_e))/abs(y_hat_rao)
    cv_dyn <- sqrt(diag(vcov_e))/abs(y_hat_dyn)
    return(list(y_rao = y_hat_rao, y_dyn = y_hat_dyn, X = X, vcov_e = vcov_e, cv_rao = cv_rao, cv_dyn = cv_dyn))
}

## Gamma_nu and Gamma_u - functions to compute gamma_nu, gamma_u, and their first derivates for
## the dynamic and Rao-Yu model
Gamma_nu <- function(rho, T, Stationarity = "FALSE") {
    Gamma_nu <- Gamma_nu.rho <- matrix(0, nrow = T, ncol = T)
    if (rho == 0) 
        return(list(Gamma_nu = Gamma_nu, Gamma_nu.rho = Gamma_nu.rho))
    if (toupper(Stationarity) == "FALSE") {
        diag(Gamma_nu) <- cumsum(rho^(2 * (0:(T - 1))))
        diag(Gamma_nu.rho) <- cumsum((2 * (0:(T - 1))) * rho^((2 * (0:(T - 1)) - 1)))
        for (j in 2:T) {
            for (i in 1:j - 1) {
                Gamma_nu[i, j] <- (rho^(j - i)) * sum(rho^(2 * (0:(i - 1))))
                Gamma_nu[j, i] <- Gamma_nu[i, j]
                Gamma_nu.rho[i, j] <- (rho^(j - i - 1)) * sum((j - i + 2 * (0:(i - 1))) * rho^(2 * 
                  (0:(i - 1))))
                Gamma_nu.rho[j, i] <- Gamma_nu.rho[i, j]
            }
        }
    } else if (toupper(Stationarity) == "TRUE") {
        diag(Gamma_nu) <- rep(1/(1 - rho^2), T)
        diag(Gamma_nu.rho) <- rep(2 * rho/(1 - rho^2)^2, T)
        for (j in 2:T) {
            for (i in 1:j - 1) {
                Gamma_nu[i, j] <- rho^(j - i)/(1 - rho^2)
                Gamma_nu[j, i] <- Gamma_nu[i, j]
                Gamma_nu.rho[i, j] <- (rho^(j - i - 1)) * (j - i + 2 * rho^2/(1 - rho^2))/(1 - rho^2)
                Gamma_nu.rho[j, i] <- Gamma_nu.rho[i, j]
            }
        }
    } else {
        print("Error: must indicate stationarity status (TRUE or FALSE)")
        break
    }
    return(list(Gamma_nu = Gamma_nu, Gamma_nu.rho = Gamma_nu.rho))
}

Gamma_u_dyn <- function(rho, T) {
    Gamma_u <- Gamma_u.rho <- matrix(0, nrow = T, ncol = T)
    if (rho == 0) 
        return(list(Gamma_u = Gamma_u, Gamma_u.rho = Gamma_u.rho))
    rho_power <- rho^(0:(T - 1))
    Gamma_u <- rho_power %*% t(rho_power)
    if (rho == 0) {
        rho_d <- rep(0, T)
    } else {
        rho_d <- c(0, (rho^(0:(T - 2))) * (1:(T - 1)))
    }
    Gamma_u.rho <- rho_power %*% t(rho_d) + rho_d %*% t(rho_power)
    return(list(Gamma_u = Gamma_u, Gamma_u.rho = Gamma_u.rho))
}

## Loglikehood function
logLik <- function(par, Model, Est.Method = "reml", y, X, vcov_e, M, T, Stationarity, randomwalk, 
    toll = 1e-122) {
    
    if (toupper(randomwalk) == "TRUE" & Stationarity == "TRUE") {
        print("Error: must model must be either stationary or random walk, not both")
        break
    }
    
    if (toupper(randomwalk) == "TRUE") {
        sig2_nu <- as.numeric(par[1])
        sig2_u <- as.numeric(par[2])
        rho <- 1
    } else if (toupper(randomwalk) == "FALSE") {
        sig2_nu <- as.numeric(par[1])
        sig2_u <- as.numeric(par[2])
        rho <- as.numeric(par[3])
    } else {
        print("Error: must indicate if random walk model (TRUE or FALSE)")
        break
    }
    
    if (sig2_nu <= 0 || sig2_u <= 0) 
        return(NA)
    
    Gamma.list <- Gamma_nu(rho = rho, T = T, Stationarity = Stationarity)
    gamma_nu <- diag(M) %x% Gamma.list$Gamma_nu
    gamma_nu.rho <- diag(M) %x% Gamma.list$Gamma_nu.rho
    
    if (tolower(Model) == "rao") {
        gamma_u <- diag(M) %x% (rep(1, T) %*% t(rep(1, T)))
        gamma_u.rho <- 0 * diag(M * T)
    } else if (tolower(Model) == "dyn") {
        Gamma.list <- Gamma_u_dyn(rho = rho, T = T)
        gamma_u <- diag(M) %x% Gamma.list$Gamma_u
        gamma_u.rho <- diag(M) %x% Gamma.list$Gamma_u.rho
    } else {
        print("Error: must indicate appropriate model (rao or dyn)")
        break
    }
    
    V <- vcov_e + sig2_nu * gamma_nu + sig2_u * gamma_u
    Xt.V.inv.X <- 0 * diag(ncol(X))
    Xt.V.inv.y <- rep(0, ncol(X))
    V.inv <- 0 * diag(M * T)
    for (d in 1:M) {
        yd <- y[((d - 1) * T + 1):(d * T)]
        Xd <- X[((d - 1) * T + 1):(d * T), ]
        Vd <- V[((d - 1) * T + 1):(d * T), ((d - 1) * T + 1):(d * T)]
        Vd.inv <- solve(Vd, tol = toll)
        Xt.V.inv.X <- Xt.V.inv.X + t(Xd) %*% Vd.inv %*% Xd
        Xt.V.inv.y <- Xt.V.inv.y + t(Xd) %*% Vd.inv %*% yd
        V.inv[((d - 1) * T + 1):(d * T), ((d - 1) * T + 1):(d * T)] <- Vd.inv
    }
    
    if (tolower(Est.Method) == "mle") {
        B_tilde <- solve(Xt.V.inv.X, Xt.V.inv.y, tol = toll)  # (6.2.5)  
        res <- y - X %*% B_tilde
        llike <- -0.5 * (determinant(V)$modulus[1] + t(res) %*% V.inv %*% res)
    } else if (tolower(Est.Method) == "reml") {
        P <- V.inv - V.inv %*% X %*% solve(Xt.V.inv.X, tol = toll) %*% t(X) %*% V.inv  # (p. 101)
        A <- diag(M * T) - X %*% solve(t(X) %*% X) %*% t(X)
        V.A <- t(A) %*% V %*% A
        # V.A.inv.P <- A%*%solve(V.A,tol=toll)%*%t(A) llike <- -0.5*(determinant(V.A)$modulus[1] +
        # t(y)%*%solve(V.A,tol=toll)%*%res)
        llike <- -0.5 * (determinant(V.A)$modulus[1] + t(y) %*% P %*% y)
    } else {
        print("Error: must indicate appropriate estimation method (reml or mle)")
        break
    }
    
    return(as.numeric(llike))
}

## Gradient of loglikehood
grad.logLik <- function(par, Model, Est.Method = "reml", y, X, vcov_e, M, T, Stationarity, randomwalk, 
    toll = 1e-122) {
    
    if (toupper(randomwalk) == "TRUE" & Stationarity == "TRUE") {
        print("Error: must model must be either stationary or random walk, not both")
        break
    }
    
    if (toupper(randomwalk) == "TRUE") {
        sig2_nu <- as.numeric(par[1])
        sig2_u <- as.numeric(par[2])
        rho <- 1
    } else if (toupper(randomwalk) == "FALSE") {
        sig2_nu <- as.numeric(par[1])
        sig2_u <- as.numeric(par[2])
        rho <- as.numeric(par[3])
    } else {
        print("Error: must indicate if random walk model (TRUE or FALSE)")
        break
    }
    
    Gamma.list <- Gamma_nu(rho = rho, T = T, Stationarity = Stationarity)
    gamma_nu <- diag(M) %x% Gamma.list$Gamma_nu
    gamma_nu.rho <- diag(M) %x% Gamma.list$Gamma_nu.rho
    
    if (tolower(Model) == "rao") {
        gamma_u <- diag(M) %x% (rep(1, T) %*% t(rep(1, T)))
        gamma_u.rho <- 0 * diag(M * T)
    } else if (tolower(Model) == "dyn") {
        Gamma.list <- Gamma_u_rho(rho = rho, T = T)
        gamma_u <- diag(M) %x% Gamma.list$Gamma_u
        gamma_u.rho <- diag(M) %x% Gamma.list$Gamma_u.rho
    } else {
        print("Error: must indicate appropriate model (rao or dyn)")
        break
    }
    
    V <- vcov_e + sig2_nu * gamma_nu + sig2_u * gamma_u
    # B_Est <- solve(Xt.V.inv.X, t(X)%*%solve(V,y,tol=toll),tol=toll) # (6.2.5)
    Xt.V.inv.X <- 0 * diag(ncol(X))
    Xt.V.inv.y <- rep(0, ncol(X))
    V.inv <- 0 * diag(M * T)
    for (d in 1:M) {
        yd <- y[((d - 1) * T + 1):(d * T)]
        Xd <- X[((d - 1) * T + 1):(d * T), ]
        Vd <- V[((d - 1) * T + 1):(d * T), ((d - 1) * T + 1):(d * T)]
        Vd.inv <- solve(Vd, tol = toll)
        Xt.V.inv.X <- Xt.V.inv.X + t(Xd) %*% Vd.inv %*% Xd
        Xt.V.inv.y <- Xt.V.inv.y + t(Xd) %*% Vd.inv %*% yd
        V.inv[((d - 1) * T + 1):(d * T), ((d - 1) * T + 1):(d * T)] <- Vd.inv
    }
    
    if (tolower(Model) == "rao") {
        V.rho <- (sig2_nu * gamma_nu.rho)
    } else if (tolower(Model) == "dyn") {
        V.rho <- (sig2_nu * gamma_nu.rho + sig2_u * gamma_u.rho)
    } else {
        print("Error: must indicate appropriate model (rao or dyn)")
        break
    }
    
    grad.llike <- rep(0, times = 3)
    
    if (tolower(Est.Method) == "mle") {
        B_tilde <- solve(Xt.V.inv.X, Xt.V.inv.y, tol = toll)  # (6.2.5)  
        res <- y - X %*% B_tilde
        V.inv.res <- V.inv %*% res
        grad.llike[1] <- -0.5 * sum(diag(V.inv %*% gamma_nu)) + 0.5 * t(V.inv.res) %*% gamma_nu %*% 
            V.inv.res
        grad.llike[2] <- -0.5 * sum(diag(V.inv %*% gamma_u)) + 0.5 * t(V.inv.res) %*% gamma_u %*% 
            V.inv.res
        grad.llike[3] <- -0.5 * sum(diag(V.inv %*% V.rho)) + 0.5 * t(V.inv.res) %*% V.rho %*% V.inv.res
    } else if (tolower(Est.Method) == "reml") {
        P <- V.inv - V.inv %*% X %*% solve(Xt.V.inv.X, tol = toll) %*% t(X) %*% V.inv  # (p. 101)
        grad.llike[1] <- -0.5 * sum(diag(P %*% gamma_nu)) + 0.5 * t(P %*% y) %*% gamma_nu %*% (P %*% 
            y)
        grad.llike[2] <- -0.5 * sum(diag(P %*% gamma_u)) + 0.5 * t(P %*% y) %*% gamma_u %*% (P %*% 
            y)
        grad.llike[3] <- -0.5 * sum(diag(P %*% V.rho)) + 0.5 * t(P %*% y) %*% V.rho %*% (P %*% y)
    } else {
        print("Error: must indicate appropriate estimation method (reml or mle)")
        break
    }
    
    if (toupper(randomwalk) == "TRUE") {
        return(as.vector(grad.llike[1:2]))
    } else if (toupper(randomwalk) == "FALSE") {
        return(as.vector(grad.llike))
    }
}

## Information matrix
info.matrix <- function(par, Model, Est.Method = "reml", y, X, vcov_e, M, T, Stationarity, randomwalk, 
    toll = 1e-22) {
    
    if (toupper(randomwalk) == "TRUE" & Stationarity == "TRUE") {
        print("Error: must model must be either stationary or random walk, not both")
        break
    }
    
    if (toupper(randomwalk) == "TRUE") {
        sig2_nu <- as.numeric(par[1])
        sig2_u <- as.numeric(par[2])
        rho <- 1
    } else if (toupper(randomwalk) == "FALSE") {
        sig2_nu <- as.numeric(par[1])
        sig2_u <- as.numeric(par[2])
        rho <- as.numeric(par[3])
    } else {
        print("Error: must indicate if random walk model (TRUE or FALSE)")
        break
    }
    
    Gamma.list <- Gamma_nu(rho = rho, T = T, Stationarity = Stationarity)
    gamma_nu <- diag(M) %x% Gamma.list$Gamma_nu
    gamma_nu.rho <- diag(M) %x% Gamma.list$Gamma_nu.rho
    
    if (tolower(Model) == "rao") {
        gamma_u <- diag(M) %x% (rep(1, T) %*% t(rep(1, T)))
        gamma_u.rho <- 0 * diag(M * T)
    } else if (tolower(Model) == "dyn") {
        Gamma.list <- Gamma_u_rho(rho = rho, T = T)
        gamma_u <- diag(M) %x% Gamma.list$Gamma_u
        gamma_u.rho <- diag(M) %x% Gamma.list$Gamma_u.rho
    } else {
        print("Error: must indicate appropriate model (rao or dyn)")
        break
    }
    
    V <- vcov_e + sig2_nu * gamma_nu + sig2_u * gamma_u
    # B_Est <- solve(Xt.V.inv.X, t(X)%*%solve(V,y,tol=toll),tol=toll) # (6.2.5)
    Xt.V.inv.X <- 0 * diag(ncol(X))
    Xt.V.inv.y <- rep(0, ncol(X))
    V.inv <- 0 * diag(M * T)
    for (d in 1:M) {
        yd <- y[((d - 1) * T + 1):(d * T)]
        Xd <- X[((d - 1) * T + 1):(d * T), ]
        Vd <- V[((d - 1) * T + 1):(d * T), ((d - 1) * T + 1):(d * T)]
        Vd.inv <- solve(Vd, tol = toll)
        Xt.V.inv.X <- Xt.V.inv.X + t(Xd) %*% Vd.inv %*% Xd
        Xt.V.inv.y <- Xt.V.inv.y + t(Xd) %*% Vd.inv %*% yd
        V.inv[((d - 1) * T + 1):(d * T), ((d - 1) * T + 1):(d * T)] <- Vd.inv
    }
    
    if (tolower(Model) == "rao") {
        V.rho <- (sig2_nu * gamma_nu.rho)
    } else if (tolower(Model) == "dyn") {
        V.rho <- (sig2_nu * gamma_nu.rho + sig2_u * gamma_u.rho)
    } else {
        print("Error: must indicate appropriate model (rao or dyn)")
        break
    }
    
    Inf_mat <- matrix(0, nrow = 3, ncol = 3)
    
    if (tolower(Est.Method) == "mle") {
        V.j.1.m <- V.inv %*% gamma_nu
        V.j.2.m <- V.inv %*% gamma_u
        V.j.3.m <- V.inv %*% V.rho
    } else if (tolower(Est.Method) == "reml") {
        P <- V.inv - V.inv %*% X %*% solve(Xt.V.inv.X, tol = toll) %*% t(X) %*% V.inv  # (p. 101)
        V.j.1.m <- P %*% gamma_nu
        V.j.2.m <- P %*% gamma_u
        V.j.3.m <- P %*% V.rho
    } else {
        print("Error: must indicate appropriate estimation method (reml or mle)")
        break
    }
    
    Inf_mat[1, 1] <- 0.5 * sum(diag(V.j.1.m %*% V.j.1.m))  # (6.2.19)
    Inf_mat[1, 2] <- 0.5 * sum(diag(V.j.1.m %*% V.j.2.m))
    Inf_mat[2, 2] <- 0.5 * sum(diag(V.j.2.m %*% V.j.2.m))
    Inf_mat[1, 3] <- 0.5 * sum(diag(V.j.1.m %*% V.j.3.m))
    Inf_mat[2, 3] <- 0.5 * sum(diag(V.j.2.m %*% V.j.3.m))
    Inf_mat[3, 3] <- 0.5 * sum(diag(V.j.3.m %*% V.j.3.m))
    Inf_mat[2, 1] <- Inf_mat[1, 2]
    Inf_mat[3, 1] <- Inf_mat[1, 3]
    Inf_mat[3, 2] <- Inf_mat[2, 3]
    
    if (toupper(randomwalk) == "TRUE") {
        return(-Inf_mat[1:2, 1:2])
    } else if (toupper(randomwalk) == "FALSE") {
        return(-Inf_mat)
    }
    
}



# reml.dyn - restricted maximum likelihood estimation of the dynamic model requires DynGamma.R


mse.est <- function(par, Model, y, X, vcov_e, M, T, Stationarity, randomwalk, Est.Method = "mle", 
    toll = 1e-22) {
    
    if (toupper(randomwalk) == "TRUE" & Stationarity == "TRUE") {
        print("Error: must model must be either stationary or random walk, not both")
        break
    }
    
    if (toupper(randomwalk) == "TRUE") {
        sig2_nu <- as.numeric(par[1])
        sig2_u <- as.numeric(par[2])
        rho <- 1
    } else if (toupper(randomwalk) == "FALSE") {
        sig2_nu <- as.numeric(par[1])
        sig2_u <- as.numeric(par[2])
        rho <- as.numeric(par[3])
    } else {
        print("Error: must indicate if random walk model (TRUE or FALSE)")
        break
    }
    
    Gamma.list <- Gamma_nu(rho = rho, T = T, Stationarity = Stationarity)
    gamma_nu <- diag(M) %x% Gamma.list$Gamma_nu
    gamma_nu.rho <- diag(M) %x% Gamma.list$Gamma_nu.rho
    
    if (tolower(Model) == "rao") {
        gamma_u <- diag(M) %x% (rep(1, T) %*% t(rep(1, T)))
        gamma_u.rho <- 0 * diag(M * T)
    } else if (tolower(Model) == "dyn") {
        Gamma.list <- Gamma_u_rho(rho = rho, T = T)
        gamma_u <- diag(M) %x% Gamma.list$Gamma_u
        gamma_u.rho <- diag(M) %x% Gamma.list$Gamma_u.rho
    } else {
        print("Error: must indicate appropriate model (rao or dyn)")
        break
    }
    
    V <- vcov_e + sig2_nu * gamma_nu + sig2_u * gamma_u
    Xt.V.inv.X <- 0 * diag(ncol(X))
    Xt.V.inv.y <- rep(0, ncol(X))
    V.inv <- 0 * diag(M * T)
    for (d in 1:M) {
        yd <- y[((d - 1) * T + 1):(d * T)]
        Xd <- X[((d - 1) * T + 1):(d * T), ]
        Vd <- V[((d - 1) * T + 1):(d * T), ((d - 1) * T + 1):(d * T)]
        Vd.inv <- solve(Vd, tol = toll)
        Xt.V.inv.X <- Xt.V.inv.X + t(Xd) %*% Vd.inv %*% Xd
        Xt.V.inv.y <- Xt.V.inv.y + t(Xd) %*% Vd.inv %*% yd
        V.inv[((d - 1) * T + 1):(d * T), ((d - 1) * T + 1):(d * T)] <- Vd.inv
    }
    
    if (tolower(Model) == "rao") {
        V.rho <- (sig2_nu * gamma_nu.rho)
    } else if (tolower(Model) == "dyn") {
        V.rho <- (sig2_nu * gamma_nu.rho + sig2_u * gamma_u.rho)
    } else {
        print("Error: must indicate appropriate model (rao or dyn)")
        break
    }
    
    Inf_mat <- matrix(0, nrow = 3, ncol = 3)
    
    if (tolower(Est.Method) == "mle") {
        V.j.1.m <- V.inv %*% gamma_nu
        V.j.2.m <- V.inv %*% gamma_u
        V.j.3.m <- V.inv %*% V.rho
    } else if (tolower(Est.Method) == "reml") {
        P <- V.inv - V.inv %*% X %*% solve(Xt.V.inv.X, tol = toll) %*% t(X) %*% V.inv  # (p. 101)
        V.j.1.m <- P %*% gamma_nu
        V.j.2.m <- P %*% gamma_u
        V.j.3.m <- P %*% V.rho
    } else {
        print("Error: must indicate appropriate estimation method (reml or mle)")
        break
    }
    
    Inf_mat[1, 1] <- 0.5 * sum(diag(V.j.1.m %*% V.j.1.m))  # (6.2.19)
    Inf_mat[1, 2] <- 0.5 * sum(diag(V.j.1.m %*% V.j.2.m))
    Inf_mat[2, 2] <- 0.5 * sum(diag(V.j.2.m %*% V.j.2.m))
    Inf_mat[1, 3] <- 0.5 * sum(diag(V.j.1.m %*% V.j.3.m))
    Inf_mat[2, 3] <- 0.5 * sum(diag(V.j.2.m %*% V.j.3.m))
    Inf_mat[3, 3] <- 0.5 * sum(diag(V.j.3.m %*% V.j.3.m))
    Inf_mat[2, 1] <- Inf_mat[1, 2]
    Inf_mat[3, 1] <- Inf_mat[1, 3]
    Inf_mat[3, 2] <- Inf_mat[2, 3]
    
    if (toupper(randomwalk) == "TRUE") {
        Inf_mat <- Inf_mat[1:2, 1:2]
    }
    
    B_Est <- solve(Xt.V.inv.X, t(X) %*% solve(V, y, tol = toll), tol = toll)  # (6.2.5)  
    
    wt1 <- wt2 <- wt3 <- rep(0, M)
    Est_fixed <- as.vector(X[T * (1:M), ] %*% B_Est)
    Est <- Est_fixed
    Xt.V.inv.X.inv <- solve(Xt.V.inv.X, tol = toll)
    M_term <- sig2_nu * gamma_nu[T, 1:T] + sig2_u * gamma_u[T, 1:T]
    for (d in 1:M) {
        r1 <- (d - 1) * T + 1
        r2 <- d * T
        Est[d] <- Est[d] + M_term %*% solve(V[r1:r2, r1:r2], (y[r1:r2] - X[r1:r2, ] %*% B_Est), tol = toll)
        tw <- M_term %*% solve(V[r1:r2, r1:r2], tol = toll)
        wt1[d] <- tw[T]
        wt2[d] <- 1 - sum(tw)
        wt3[d] <- sum(tw[1:(T - 1)])
    }
    
    # begin MSE calculation
    g1 <- g2 <- g3 <- rep(0, M)
    mi <- c(rep(0, T - 1), 1)
    b <- matrix(0, nrow = T, ncol = M)
    Gi <- sig2_nu * gamma_nu[1:T, 1:T] + sig2_u * gamma_u[1:T, 1:T]
    G <- diag(M) %x% (Gi)
    V_bar <- solve(Inf_mat, tol = toll)
    for (d in 1:M) {
        Xi <- X[((d - 1) * T + 1):(d * T), ]
        li <- as.vector(Xi[T, ])
        Vi <- V[((d - 1) * T + 1):(d * T), ((d - 1) * T + 1):(d * T)]
        vcov_ei <- vcov_e[((d - 1) * T + 1):(d * T), ((d - 1) * T + 1):(d * T)]
        Vi.inv <- solve(Vi, tol = toll)
        bi <- as.vector(t(mi) %*% Gi %*% Vi.inv)  # (6.2.10)
        di <- li - as.vector(t(bi) %*% Xi)  # (6.2.10)
        if (toupper(randomwalk) == "FALSE") {
            d.bi <- matrix(0, nrow = T, ncol = 3)
        } else {
            d.bi <- matrix(0, nrow = T, ncol = 2)
        }
        d.bi[, 1] <- t(t(mi) %*% gamma_nu[1:T, 1:T] %*% Vi.inv - t(mi) %*% Gi %*% Vi.inv %*% gamma_nu[1:T, 
            1:T] %*% Vi.inv)
        d.bi[, 2] <- t(t(mi) %*% gamma_u[1:T, 1:T] %*% Vi.inv - t(mi) %*% Gi %*% Vi.inv %*% gamma_u[1:T, 
            1:T] %*% Vi.inv)
        if (toupper(randomwalk) == "FALSE") {
            d.Gi_rho <- (sig2_nu * gamma_nu.rho[1:T, 1:T] + sig2_u * gamma_nu.rho[1:T, 1:T])
            d.bi[, 3] <- t(t(mi) %*% d.Gi_rho %*% Vi.inv - t(mi) %*% Gi %*% Vi.inv %*% (sig2_nu * 
                gamma_nu.rho[1:T, 1:T] + sig2_u * gamma_nu.rho[1:T, 1:T]) %*% Vi.inv)
        }
        g1[d] <- t(mi) %*% (Gi - Gi %*% Vi.inv %*% Gi) %*% mi
        g2[d] <- t(di) %*% Xt.V.inv.X.inv %*% di
        g3[d] <- sum(diag((t(d.bi) %*% Vi %*% d.bi) %*% V_bar))
    }
    
    wt <- data.frame(wt1, wt2, wt3)
    names(wt) <- c("resid.weight", "synthetic.weight", "past.resid.weight")
    mse <- g1 + g2 + 2 * g3
    mse.decomp <- data.frame(g1, g2, g3)
    return(list(coef = B_Est, Est = Est, Est_fixed = Est_fixed, compositing.weights = wt, mse = mse, 
        mse.decomp = mse.decomp))
}
