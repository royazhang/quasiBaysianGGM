####### Graphical Bayesian Lasso #######

library(MASS)
library(GIGrvg)

################################
BGL <- function(Y, chain_len, burnin) {

  n = dim(Y)[1]
  p = dim(Y)[2]
  S = t(Y)%*%Y
  
  W_path = array(data=NA, dim=c(p,p,(chain_len-burnin)/10))
  
  # initial values
  W = diag(5,p)
  lambda = 10 # parameter lambda ~ Gamma(shape=1,rate=0.01)
  tau = matrix(1, p, p)
  
  ##### Gibbs Sampler #####
  for (iter in 1:chain_len) {
    for (i in 1:p) {
      # sample ith column in Omega
      invOmega11 =  chol2inv(chol(W[-i,-i]))
      r = rgamma(1, shape=n/2+1, rate=(S[i,i]+lambda)/2)
      C = chol2inv(chol((S[i,i] + lambda)*invOmega11 + diag(1/tau[i,-i])))
      #Omega11_s = W[-i,-i]/(S[i,i]+lambda)
      #C = Omega11_s %*% (diag(1,p-1)-solve(diag(tau[i,-i]) + Omega11_s) %*% Omega11_s)
      beta = mvrnorm(1, mu=-C%*%S[i,-i], Sigma=C) 
      W[i,-i] = beta
      W[-i,i] = beta
      W[i,i] = r + beta %*% invOmega11 %*% beta
      
      # sample tau
      tau[i,-i] = 1/mapply(rgig, n=1, lambda=rep(-0.5,p-1), psi=beta^2, chi=rep(lambda^2,p-1))
      tau[-i,i] = tau[i,-i]
    }
    # sample lambda
    lambda = rgamma(1, shape=1+p*(p+1)/2, rate=0.01+sum(abs(W))/2)
    
    if (iter>burnin & iter%%10 == 0) {
      W_path[,,(iter-burnin)/10] = W
    }
  }
  #####
  return(W_path)
}
################################



