# functions used in Step 2 for updating z

zbar_gaussian <- function(
  X, # n by np
  Y, # n by 1
  beta.m1, # np by 1
  theta.m, # n by 1
  sigma2, # variance, probability need to estimate
  rho
  ){
  
  n <- dim(X)[1]
  p <- dim(X)[2]/n
  
  a_bar <- X %*% beta.m1/p
  z_bar.m1 <- (Y + sigma2*rho * (a_bar + theta.m))/(p + sigma2*rho)
  return(z_bar.m1)
}

# zbar_gaussian(X, Y, beta.m1 = beta, theta.m = rep(1, 200),
#               sigma2 = 1, rho = .1)
