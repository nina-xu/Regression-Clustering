# 3/7/2019
### Step 3: update theta
update_theta <- function(
  X,
  beta.m1, # from this iteration
  theta.m, # from previous iteration
  z_bar.m1 # from this iteration
){
  n <- dim(X)[1]
  p <- dim(X)[2]/n
  theta.m1 <- theta.m + X %*% beta.m1/p - z_bar.m1 
  #output
  return(theta.m1)
}

# # test
# z_bar.m1 <- zbar_gaussian(X, beta.m1 = beta, theta.m = rep(1, 200),
#               sigma2 = 1, rho = .1)
# update_theta(X, Y, beta.m1 = beta, theta.m = rep(1, 200), z_bar.m1 = z_bar.m1)
