library(Matrix)
library(ggplot2)
library(scales)
library(gridExtra)
library(GPDAG)

## model and priors
nu = 1.5
beta = 1.5
sig = 0.1
sig_prior_size = 10
nug_prior = c(sig_prior_size, sig_prior_size*sig^2)   ## gamma distribution parameters for sig2^{-1}
tau = 10
tau_bound = c(tau-1e-3,tau+1e-3)
log_tau_fun <- function(tau){
  return(1)
}

## dataset; build norming dag
set.seed(1110)
J = 8
n = 2^J+1

X = seq(0,1,1/(n-1))
Y = runif(n)
X_cov = cov_matern(matrix(X), matrix(X), nu, tau)
Yf = MASS::mvrnorm(1,rep(0,n),X_cov)
Y = Yf + rnorm(n)*sig

sig2_bound = c(1e-3, sqrt(sum((Y-mean(Y))^2)/(n-1)) )

## Experimental settings
n_mcmc = 2000
n_burn = 1000
j_lst = seq(4,J)
nj = length(j_lst)
fhat_records = vector('list',nj)
finf_records = rep(0,nj)
sig_records = rep(0,nj)
fhat_mm_records = vector('list',nj)
finf_mm_records = rep(0,nj)
sig_mm_records = rep(0,nj)
df_records = vector('list',nj)
plots = vector('list',nj)
RunTime_records = matrix(0,nrow=nj,ncol=2)
for (j in j_lst){
  j_ids = seq(1,n,2^(J-j))
  mcmc_obj = GPgrid(Y[j_ids],minsep=1/2^j, adapt=FALSE, nu=nu, sig=sig,
                    tau=tau, nug_prior=nug_prior, n_mcmc=2000, n_burn=1000)
  fhat_obj = confidence_bands(mcmc_obj$Z_mcmc)
  sighat = mean(mcmc_obj$sig_mcmc)
  print(paste('n=',2^j+1))
  print(paste('l^infty estimation error, Norming DAG:', max(abs(fhat_obj$mean-Yf[j_ids][mcmc_obj$sort_in_ord]))))

  ## records mcmc outputs
  j_ind = j-j_lst[1]+1
  fhat_records[[j_ind]] = fhat_obj$mean
  finf_records[j_ind] = max(abs(fhat_obj$mean-Yf[1:(2^j+1)]))
  sig_records[j_ind] = sighat
}

