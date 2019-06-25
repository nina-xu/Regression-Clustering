
objective_gaussian <- function(X, Y, beta, lambda, nu, B,
                                weight.w, # as a list
                                weight.u, sigma2,
                                variable.selection = F){
  n <- length(Y); 
  p <- ncol(X)/n
  obj1 <- 1/2/sigma2*sum((Y - X %*% beta)^2)
  beta.matrix <- matrix(beta, n, p)
  Bb <- as.vector(B %*% beta.matrix)
  obj2 <- lambda * t(unlist(weight.w)) %*% abs(Bb) 
  obj <- obj1 + obj2
  # need to add penalty term for variable selection later
  if(variable.selection == T){
    weightu.matrix <- matrix(unlist(weight.u), n, p)
    obj3 <- nu * t(unlist(weight.u)) %*% as.vector(abs(beta)) 
    obj <- obj +  obj3
    return(list(obj=obj, obj1.likeli = obj1, obj2.pen.lambda = obj2, obj3.pen.nu = obj3))
  } else{
    return(list(obj=obj, obj1.likeli = obj1, obj2.pen.lambda = obj2))
  }
}
# objective_gaussian <- function(X, # n by np
#                                Y, # n by 1
#                                beta, # np by 1
#                                lambda, nu, 
#                                weight.w, # as a list
#                                weight.u, # as a list
#                                sigma2,
#                                #B= B,
#                                variable.selection = F){
#   # n <- length(Y); p <- ncol(X)/n
#   eks = ek(n)
#   # B = t( eks$ek1- eks$ek2) # n2-by-n dim
#   l <- -1/2/sigma2*sum((Y - X %*% beta)^2)
#   obj1 <- -l/n
#   beta.matrix <- matrix(beta, n, p)
#   Bb <- as.vector(B %*% beta.matrix)
#   obj2 <- t(unlist(weight.w)) %*% abs(Bb) 
#   obj <- obj1 + lambda * obj2
#     # need to add penalty term for variable selection later
#   if(variable.selection == T){
#     weightu.matrix <- matrix(unlist(weight.u), n, p)
#     obj3 <- t(unlist(weight.u)) %*% as.vector(abs(beta)) 
#     obj <- obj + nu * obj3
#     return(list(obj=obj, obj1.likeli = obj1, obj2.pen.lambda = obj2, obj3.pen.nu = obj3))
#   } else{
#     return(list(obj=obj, obj1.likeli = obj1, obj2.pen.lambda = obj2))
#     }
# }


