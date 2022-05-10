# gLasso
library(CVglasso)
library(MASS)

setwd("C:\\Users\\zhang\\documents\\new_sim\\run_hub_clique")
load("Omega0_hubclq_p100n50.Rdata")
#Omega0 = Omega0_clq

p = 100
n = 150 
nDataSet = 50
gLasso_pstmean50 = array(NA, dim=c(p,p,nDataSet))
gLasso_isZero50 = array(NA, dim=c(p,p,nDataSet))
gLasso_cputime50 = rep(0, nDataSet)


for (i in 1:nDataSet) {
  Y = mvrnorm(n, mu=rep(0,p), Sigma=chol2inv(chol(Omega0)))
  
  pt = proc.time()
  result = CVglasso(X=Y, diagonal=FALSE)
  alltime = proc.time()-pt
  gLasso_pstmean50[,,i] = result$Omega
  gLasso_isZero50[,,i] = (result$Omega==0)
  gLasso_cputime50[i] = alltime[1]
}

# summary result
Fnorm_mean = apply(gLasso_pstmean50, 3, function(x) sqrt(sum((x-Omega0)^2)))
mean(Fnorm_mean)
sd(Fnorm_mean)

getTPR = function(x) (sum(Omega0!=0 & !x)-100)/(sum(Omega0!=0)-100)
getFPR = function(x) sum(Omega0==0 & !x)/sum(Omega0==0)

TPR = apply(gLasso_isZero50, 3, getTPR)
mean(TPR)*100
sd(TPR)*100

FPR = apply(gLasso_isZero50, 3, getFPR)
mean(FPR)*100
sd(FPR)*100

mean(gLasso_cputime50)/20
