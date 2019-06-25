
### main function
ADMM <- function(
  X0, # n * p
  Y, # n * 1
  lambda, # tuning parameter for fused-lasso penalty (clustering)
  nu=0, # tuning parameter for lasso penalty (variable selection)
  rho, # ascent rate for dual variable theta, I believe
  family = 'gaussian',
  variable.selection = F, # right now we are not doing variable selection
  maxite = 1000, # maximum number of iterations
  epri = 0.1, # stopping criterion for primal residual
  edual = 0.1, # stopping criterion for dual residual
  weight.w = split(rep(1, p*n*(n-1)/2), rep(1:p, each=n*(n-1)/2)), # weights for fusing
  weight.u = split(rep(1, n*p), rep(1:p, each=n)), # weights for variable selection
  packages = "doParallel",
  export = c("general.fusedlasso.solver","general.fusedlasso.lasso.solver",
             "proximal.L1")
){
  n <- nrow(X0)
  p <- ncol(X0)
  # rescale weight
  weight.w = lapply(weight.w,function(data) data/sum(abs(data)))
  weight.u = lapply(weight.u,function(data) data/sum(abs(data)))
  
  if (n != length(Y)) stop("dimensions of X and Y don't match")
  converge = 0
  
  # Initialize
  initial.values <- initial(X0, Y, family = family)
  
  X <- initial.values$X.matrix
  X.list <- initial.values$X.list
  Y <- initial.values$Y
  B <- initial.values$B
  beta <- beta0 # this is just for simulation. Using true value as initial value
  # beta <- initial.values$beta.0
  z_bar0 <- initial.values$z_bar.0
  theta <- initial.values$theta.0 
  sigma2 <- initial.values$sigma2
  
  # if(family == 'gaussian'){
  #   sigma2 <- initial.values$sigma2
  # }
  
  # Store results
  Allbeta <- list()
  Allzbar <- list()
  Alltheta <- list()
  all.fused.lasso.converge = numeric()
  
  rl2 <- numeric() # L2 norm of primal residual
  sl2 <- numeric() # L2 norm of dual residual
  # Looping over three alternating steps
  # t0 = Sys.time()
  for(m in 1:maxite){ #stopping criterion? See admm tutorial section 3.3
    # print(m)
    # Parallel step 1
    
    step1results <- update_beta(z.m = z_bar0, theta.m = theta, 
                                beta.m = beta,
                                X.list = X.list,
                                weight.w = weight.w,
                                weight.u = weight.u,
                                lambda = lambda, 
                                nu = nu,
                                rho = rho,
                                B = B,
                                variable.selection = variable.selection,
                                max.ite = 1000,
                                tol = 1e-4,
                                packages = packages,
                                export = export)
    beta = lapply(step1results,function(data) return(data$beta.m1))
    beta <- unlist(beta)
    Allbeta[[m]] <- beta
    fused.lasso.converge = lapply(step1results,function(data) return(data$converge))
    all.fused.lasso.converge[m] = mean(unlist(fused.lasso.converge))
    
    
    # Shared step 2
    # for now, only for Gaussian
    z_bar <- zbar_gaussian(X, Y, 
                           beta.m1 = beta,
                           theta.m = theta,
                           sigma2 = sigma2,
                           rho = rho
    )
    Allzbar[[m]] <- z_bar
    
    
    # update step 3
    theta <- update_theta(X, beta.m1 = beta, theta.m = theta, z_bar.m1 = z_bar)
    Alltheta[[m]] <- theta
    
    # calculate residual for stopping criteria
    rl2[m] <- sqrt(sum((X %*% beta/p - z_bar)^2/n))
    sl2[m] <- sqrt(sum((z_bar - z_bar0)^2/n))
    if(rl2[m] <= epri & sl2[m] <= edual) {
      converge <- 1
      break
    }
    
    z_bar0 <- z_bar
  }
  # Sys.time() - t0
  # closeAllConnections()
  # calculate loss
  All.loss <- objective_gaussian(X, Y, beta, lambda, nu,B,
                                 weight.w, weight.u, sigma2,variable.selection)
  
  #output
  return(list(beta = matrix(beta, n,p), Allbeta=Allbeta, Allzbar=Allzbar, Alltheta=Alltheta, Allloss = All.loss,
              fused.lasso.converge = mean(all.fused.lasso.converge),converge = converge, n.iterations = m))
}