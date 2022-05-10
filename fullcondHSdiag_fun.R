####### Graphical Quasi-Bayesian: Horseshoe #######
library(MASS)
library(glmnet)

################################
fullcondHSdiag <- function(Y, chain_len, burnin) {
  
  n = dim(Y)[1]
  p = dim(Y)[2]
  S = t(Y)%*%Y
  
  W_path = array(data=NA, dim=c(p,p,(chain_len-burnin)/10))
  
  # empirical estimates of W_ii
  sigma2 = rep(0, p)
  for (i in 1:p) {
    cv.lambda = cv.glmnet(x=Y[,-i], y=Y[,i], nfolds=10, family="gaussian", alpha=1)$lambda.min
    temp = glmnet(x=Y[,-i], y=Y[,i], family="gaussian", alpha=1, lambda=cv.lambda)
    sigma2[i] = sum((Y[,i] - Y[,-i]%*%temp$beta)^2)/(n-temp$df)
  }
  
  W = diag(1/sigma2,p)
  
  lambda = matrix(1, p, p)
  v = lambda
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
    
    temp_sum = W^2/lambda
    diag(temp_sum) = 0
    
    tau2 = 1/rgamma(n=1, shape=(p*(p-1)+1)/2, rate=1/kappa+sum(temp_sum)/2)
    kappa = 1/rexp(n=1, rate=1+1/tau2)
  }
  #####
  return(W_path)
}
################################
