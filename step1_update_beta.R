# functions used in Step 1 for updating beta

# packages = c("MASS")



# # generate B
# eks = ek(n)
# B = t( eks$ek1- eks$ek2) # n2-by-n dim

update_beta <- function(
  z.m, # from previous iteration
  theta.m, # from previous iteration
  beta.m, # from previous iteration, np by 1
  X.list,
  weight.w,
  weight.u,
  lambda,
  nu,
  rho, # constant for augmented Lagrandian function
  B, # matrix for fused items
  variable.selection = F,
  max.ite = 200,
  tol = 1e-4,
  packages = "doParallel",
  export = c("general.fusedlasso.solver","general.fusedlasso.lasso.solver",
             "proximal.L1")
){
  p = length(X.list)
  X = do.call(cbind, X.list)
  # beta.m = do.call(rbind, list.beta.m)
  list.beta.m <- split(beta.m, rep(1:p, rep(n,p)))
  bar.x.beta.m = X %*% beta.m / p
  
  # fused lasso problem
  registerDoParallel(cores=4)
  
  fused.lasso.result = foreach(j = 1:p, .packages = packages,
                               .export = export) %dopar% {

    v.j.m = as.vector(X.list[[j]] %*% list.beta.m[[j]] + z.m - bar.x.beta.m - theta.m) # pseudo outcome
    # fused.lasso.fit = genlasso(v.j.m, list.X[[j]], B) # can only handle fused lasso part with weight.w=1
    # beta.m1 = coef(fused.lasso.fit, lambda = lambda/rho)
    
    # 3/10/2019, add option of variable selection
    if(variable.selection == F)
    {
      fused.lasso.fit =general.fusedlasso.solver(Y=v.j.m,X=X.list[[j]],B,
                                                gamma1=lambda/rho,
                                                weight.w[[j]], 
                                               rho,initial.beta=list.beta.m[[j]],
                                          max.ite = max.ite,tol = tol,output = 0)
      print(fused.lasso.fit$converge)
    } else {
      fused.lasso.fit =general.fusedlasso.lasso.solver(v.j.m,X.list[[j]],B,
                                                 gamma1=lambda/rho,gamma2=nu/rho,
                                                 weight.w[[j]], 
                                                 weight.u[[j]],
                                                 rho,initial.beta=list.beta.m[[j]],
                                                 max.ite = max.ite,tol = tol,output = 0)
      print(fused.lasso.fit$converge)
    }
    list(beta.m1 = fused.lasso.fit$beta, converge = fused.lasso.fit$converge, 
         niter=fused.lasso.fit$niter)
    
  }
  list.beta.m1 = fused.lasso.result
  closeAllConnections()
  # output beta_j^m+1
  return(list.beta.m1)
}

