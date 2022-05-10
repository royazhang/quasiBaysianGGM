####### Graphical Quasi-Bayesian: Horseshoe #######
library(MASS)

################################
fullcondHS <- function(Y, chain_len, burnin) {
  
  n = dim(Y)[1]
  p = dim(Y)[2]
  S = t(Y)%*%Y
  
  W_path = array(data=NA, dim=c(p,p,(chain_len-burnin)/10))
  
  # initial values
  W = diag(1,p)
  lambda = matrix(1, p, p)
  v = matrix(1, p, p)
  tau2 = 1
  kappa = 1
  
  ##### Gibbs Sampler #####
  for (iter in 1:chain_len) {
    for (i in 1:p) {
      for (j in c(1:p)[-i]) {
        num = W[-j,i] %*% S[-j,j]
        denom = S[j,j] + W[i,i]/tau2/lambda[j,i]
        W[j,i] = rnorm(n=1, mean = -num/denom, sd = sqrt(W[i,i]/denom))
      }
    }
    
    if (iter > burnin & iter%%10 == 0) {
      W_path[,,(iter-burnin)/10] = W
    }
    
    lambda = 1/rexp(n = p^2, rate = as.vector(W^2/2/tau2 + 1/v))
    v = 1/rexp(n = p^2, rate = 1 + 1/lambda)
    lambda = matrix(lambda, ncol = p, nrow = p)
    v = matrix(v, ncol = p, nrow = p)
    
    temp_mat = W^2/lambda
    diag(temp_mat) = 0
    
    tau2 = 1/rgamma(n=1, shape=(p*(p-1)+1)/2, rate=1/kappa+sum(temp_mat)/2)
    kappa = 1/rexp(n=1, rate=1+1/tau2)
  }
  #####
  return(W_path)
}
################################
