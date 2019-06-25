#3/8/2019

# generate matrix B
ek = function(n){
  # n,k
  count = 0
  ek1 = ek2 = matrix(0,n,n*(n-1)/2)
  
  for(i in (n-1):1){
    temp = matrix(0,n,i)
    temp[n-i,] = 1
    ek1[,(count+1):(count+i)] = temp
    ek2[(n-i+1):n,(count+1):(count+i)] = diag(1,i,i)
    count = count +i
  }  
  
  return(list(ek1 = ek1, ek2 = ek2))
}

### Step 0: initialization
# input: X0, Y
# output: rearranged X, theta, z
initial <- function(X0, Y, family = 'gaussian', center = F){
  
  # centering X0 and Y
  if(center == T){
  X0 <- scale(X0, center=T)
  Y <- scale(Y, center = T)
  }
  # organize X into the correct dimension
  X <- matrix(apply(X0, 2, diag), n, n*p)
  # Makre X into a list
  X0.list <- split(X0, rep(1:p, each = n))
  X.list <- lapply(X0.list, diag)
  
  
  # Initializing beta
  # Initializing all betas to be the same, using a centered regular regression
  initialmodel <- glm(Y ~ 0 + X0, family = family)
  beta.0 <- rep(initialmodel$coefficients, rep(n, p))
  if(family == 'gaussian'){
    sigma2 <- sigma(initialmodel)^2
  }
  
  # Initializing multiplier theta
  theta.0 <- rep(1, n)  
  
  # Initializing z bar
  # Question: can initial z bar be calculated from initial beta?
  z_bar.0 <- X %*% beta.0/p
  
  # Generating B
  eks = ek(n)
  B = t( eks$ek1- eks$ek2) # n2-by-n dim
  
  if(family == 'gaussian'){
    return(list(X.matrix=X, X.list = X.list, Y=Y, B = B, beta.0 = beta.0, z_bar.0 = z_bar.0, 
                theta.0 = theta.0, sigma2 = sigma2))
  }
  
}

# # test
# X = X0
# rho = 0.1
# initialvalues <- initial(X, Y)
# dim(initialvalues$Xc)

# an example weight.w function for known cluster structure: first half belongs to one cluster 
# and the rest belongs to another. A smarter way is to make a big matrix for cluster structure and take
# the lower triangle
pseudo.weight.w = function(n,p){
  weight.w.j = numeric()
  
  for (i in 1:(n/2-1)){
    weight.w.j = c(weight.w.j, rep(1,n/2-i),rep(0,n/2))
  }
  weight.w.j = c(weight.w.j,rep(0,n/2))
  for (i in (n/2+1):(n-1)){
    weight.w.j = c(weight.w.j,rep(1,n-i))
  }
  weight.w.j = weight.w.j/sum(weight.w.j) # sum to 1
  return(split(rep(weight.w.j,p), rep(1:p, each=n*(n-1)/2)) )
  
}

pseudo.weight.w1 = function(n,p){

  cluster.structure = as.matrix(bdiag(matrix(1,n/2,n/2),matrix(1,n/2,n/2)))
  cluster.structure.lower = cluster.structure[lower.tri(cluster.structure, diag = FALSE)]
  cluster.structure.lower = as.vector(cluster.structure.lower)
  cluster.structure.lower = cluster.structure.lower/sum(cluster.structure.lower) # sum to 1
  return(split(rep(cluster.structure.lower,p), rep(1:p, each=n*(n-1)/2)) )
  
}