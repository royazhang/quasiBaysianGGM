######### run condHS #########print
source("fullcondHS_fun.R")
load("Omega0_p100_rnd.Rdata")

p = 100
n = 80
chain_len = 6000
burnin = 1000

Y = mvrnorm(n, mu=rep(0,p), Sigma=chol2inv(chol(Omega0)))
pt <- proc.time()
W_path = fullcondHS(Y, chain_len, burnin)
alltime = proc.time()-pt
condHS_rnd_cputime = alltime[1]

# point estimator
condHS_rnd_est_mean = apply(W_path, c(1,2), mean)
condHS_rnd_est_med = apply(W_path, c(1,2), median)

# selection
W_low = apply(W_path, c(1,2), quantile, 0.25, type=4)
W_high = apply(W_path, c(1,2), quantile, 0.75, type=4)
condHS_rnd_isZero = (W_low<0) & (0<W_high)

filename = paste("condHS_rnd_", Sys.getenv("SLURM_ARRAY_TASK_ID"), sep="")
save(condHS_rnd_est_mean, condHS_rnd_est_med, condHS_rnd_isZero, condHS_rnd_cputime, 
     file = paste(filename,".Rdata",sep=""))
