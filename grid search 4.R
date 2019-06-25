#3/10/2019
# grid search 2
# search for lambda 10^(-5) to 1
# using true value as initial value
setwd('G:\\My Drive\\Regression Clustering with Prof. Nina Ning Xu\\simulation')
# setwd("C:\\Users\\wangb09\\Google Drive\\Research\\Regression Clustering with Prof. Nina Ning Xu\\simulation")

#install.packages('doParallel')
library(foreach)
library(doParallel)
library(Matrix)
packages = c("MASS", "doParallel","Matrix")
export = c("general.fusedlasso.solver","general.fusedlasso.lasso.solver",
          "proximal.L1")

source("initial.R")
source("step1_update_beta.R")
source("step2_zbar_gaussian.R")
source("step3_update_theta.R")
source("gll.solver.R")
source("main function.R")
source("objective_gaussian.R")

# test
set.seed(3062019)
n <- 100
p <- 3
beta1 <- c(rep(1, n/2), rep(3, n/2))
beta2 <- c(rep(1, n/2), rep(-1, n/2))
beta3 <- c(rep(1, n/2), rep(1, n/2))
X0 <- matrix(rnorm(n*p), n, p)
e <- rnorm(n,0,1)
X <- cbind(diag(X0[,1]), diag(X0[,2]), diag(X0[,3]))
beta0 <- c(beta1, beta2, beta3)
Y <- X %*% beta0 + e
data <- cbind(Y, X0)

lambdas <- exp(seq(-7, 0, length.out = 100))

fakeweight <- rep(1,n*(n-1)/2*p)
fakeweight <- split(fakeweight, rep(1:p, rep(n*(n-1)/2,p)))
fakeweight <- lapply(fakeweight,function(data) data/sum(abs(data)))

times <- numeric()
gridsearch3betas <- list()
gridsearch3loss <- list()

for(i in 1:length(lambdas)){
  t0 <- Sys.time()
  test <- ADMM(
    X0, # n * p
    Y, # n * 1
    lambda = lambdas[i], # tuning parameter for fused-lasso penalty (clustering)
    nu = 0, # tuning parameter for lasso penalty (variable selection)
    rho = .01, # ascent rate for dual variable theta, I believe
    family = 'gaussian',
    variable.selection = F, # right now we are not doing variable selection
    maxite = 200, # maximum number of iterations
    epri = 1e-3, # stopping criterion for primal residual
    edual = 1e-3, # stopping criterion for dual residual
    weight.w = fakeweight,
    packages=packages,
    export = export
  )
  gridsearch3betas[[i]] <- test$beta
  gridsearch3loss[[i]]<- test$Allloss
  print(gridsearch3betas[[i]])
  print(test$fused.lasso.converge)
  times[i] <- Sys.time() - t0
  print(times[i] )
}

loss <- matrix(unlist(gridsearch3loss), length(lambdas), 3, byrow = T)
colnames(loss) <- c('Overall', '-loglikeli/n', 'penalty.lambda')
plot(lambdas,loss[,1], ylim = c(min(loss), max(loss)), type='l')
lines(lambdas,loss[,2],col=2)
lines(lambdas,loss[,3],col=3)
legend('topleft',colnames(loss), col=1:3, lty=1)

write.csv(cbind(lambdas, loss), 'Loss function tradeoff 031219.csv', row.names = F)

# plot for all lambdas, final results
pdf(file = 'grid search 4 all 1 weight.pdf', width = 10, height = 4)
for(i in 1:length(lambdas)){
  par(mfrow=c(1,3))
  hist(gridsearch3betas[[i]][,1], main = paste('beta1, lambda =', round(lambdas[i],5)))
  hist(gridsearch3betas[[i]][,2], main = paste('beta2, lambda =', round(lambdas[i],5)))
  hist(gridsearch3betas[[i]][,3], main = paste('beta3, lambda =', round(lambdas[i],5)))
}
dev.off()
