# gSCAD
library(GGMncv)
library(MASS)
setwd("C:\\Users\\zhang\\documents\\new_sim\\run_hub_clique")
load("Omega0_hubclq_p100n50.Rdata")
#Omega0 = Omega0_clq

p = 100
n = 150 
nDataSet = 50

gSCAD_pstmean50 = array(NA, dim=c(p,p,nDataSet))
gSCAD_isZero50 = array(NA, dim=c(p,p,nDataSet))
gSCAD_cputime50 = rep(0, nDataSet)

for (i in 1:nDataSet) {
  Y = mvrnorm(n, mu=rep(0,p), Sigma=chol2inv(chol(Omega0)))
  
  pt = proc.time()
  result = ggmncv(R=cor(Y), n, penalty="scad", initial = diag(5, p))
  alltime = proc.time()-pt
  gSCAD_pstmean50[,,i] = result$Theta
  gSCAD_isZero50[,,i] = (result$Theta==0)
  gSCAD_cputime50[i] = alltime[1]
}

# summary result
Fnorm_mean = apply(gSCAD_pstmean50, 3, function(x) sqrt(sum((x-Omega0)^2)))
mean(Fnorm_mean)
sd(Fnorm_mean)

getTPR = function(x) (sum(Omega0!=0 & !x)-100)/(sum(Omega0!=0)-100)
getFPR = function(x) sum(Omega0==0 & !x)/sum(Omega0==0)

TPR = apply(gSCAD_isZero50, 3, getTPR)
mean(TPR)*100
sd(TPR)*100

FPR = apply(gSCAD_isZero50, 3, getFPR)
mean(FPR)*100
sd(FPR)*100

mean(gSCAD_cputime50)/20
