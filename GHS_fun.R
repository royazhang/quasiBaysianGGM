####### Graphical Horseshoe #######
library(MASS)
library(GIGrvg)

################################
BGHS <- function(Y, chain_len, burnin) {
  #####
  n = dim(Y)[1]
  p = dim(Y)[2]
  S = t(Y)%*%Y
  
  W_path = array(data=NA, dim=c(p,p,(chain_len-burnin)/10))
  # initial values
  W = diag(5,p)
  lambda = matrix(1, p, p) # parameter lambda ~ Gamma(shape=1,rate=0.01)
  nu = lambda
  tau2 = 1
  eta = 1
  
  # Gibbs Sampling
  for (iter in 1:chain_len) {
    for (i in 1:p) {
      r = rgamma(n=1, shape=n/2+1, rate=S[i,i]/2)
      invOmega11 = chol2inv(chol(W[-i,-i]))
      C = chol2inv(chol(S[i,i]*invOmega11 + diag(1/(tau2*lambda[-i,i]))))
      beta = mvrnorm(1, mu=-C%*%S[-i,i], Sigma=C)
      W[-i,i] = beta
      W[i,-i] = beta
      W[i,i] = r + beta %*% invOmega11 %*% beta
      lambda[-i,i] = 1/rgamma(n=p-1, shape=rep(1,p-1), rate=1/nu[-i,i]+W[-i,i]^2/2/tau2)
      lambda[i,-i] = lambda[-i,i]
      nu[-i,i] = 1/rgamma(n=p-1, shape=rep(1,p-1), rate=1+1/lambda[-i,i])
      nu[i,-i] = nu[-i,i]
    }
    tau2 = 1/rgamma(n=1, shape=(p*(p-1)/2+1)/2, rate=1/eta+sum((W^2/lambda)[upper.tri(W, diag=FALSE)])/2)
    eta = 1/rgamma(n=1, shape=1, rate=1+1/tau2)
    
    if (iter>burnin & iter%%10 == 0) {
      W_path[,,(iter-burnin)/10] = W
    }
  }
  #####
  return(W_path)
}
################################


