library(Matrix)
library(ggplot2)
library(scales)
library(gridExtra)
library(GPDAG)
source("/home/yichen/GPDAG/paper/paper_utils.R")

##------------------------ truth is Sobolev 1/2 smooth -------------------------
## generate the underlying truth
beta = 1.5
set.seed(7)
J = 11
n = 2^J+1
X = matrix((0:(n-1)) / (n-1), nrow=n, ncol=1)
tau = 10
sig = 0.1
X_cov = cov_matern(X, X, beta, tau)
Yf = MASS::mvrnorm(1,rep(0,n),X_cov)
Yf = Yf - mean(Yf)
Y = Yf + rnorm(n)*sig

## shared prior settings for two Norming methods
sig_prior_size = 10
sig2_prior = c(sig_prior_size, sig_prior_size*sig^2)   ## gamma distribution parameters for sig2^{-1}
tau_bound_const = c(tau-1e-3,tau+1e-3)
tau_bound_adapt = c(0.1, 1000)
sig2_bound = c(1e-3, sqrt(sum((Y-mean(Y))^2)/(n-1)) )

## Norming methods. Notice different norming methods share the same layer structures, but with different DAGs
nu1 = 1.5
nu2 = 2.5
nu3 = 2.5
dag_norm1_obj = DAGgrid_per_1d(J, nu1)
X_norm = matrix(dag_norm1_obj$X_ord,nrow=n,ncol=1)
norm_ord_2_coord_ord = sort(X_norm,index.return=TRUE)$ix
coord_ord_2_norm_ord = sort(norm_ord_2_coord_ord,index.return=TRUE)$ix
Y_norm = Y[coord_ord_2_norm_ord]
Yf_norm = Yf[coord_ord_2_norm_ord]

dag_norm1 = dag_norm1_obj$dag_ord
dag_norm1_cpp = Rdag_to_Cppdag(dag_norm1)

dag_norm2_obj = DAGgrid_per_1d(J, nu2)
dag_norm2 = dag_norm2_obj$dag_ord
dag_norm2_cpp = Rdag_to_Cppdag(dag_norm2)

dag_norm3_obj = DAGgrid_per_1d(J, nu3)
dag_norm3 = dag_norm3_obj$dag_ord
dag_norm3_cpp = Rdag_to_Cppdag(dag_norm3)


## Experimental settings
n_mcmc = 4000
n_burn = 2000
j_lst = seq(4,J)
nj = length(j_lst)

fhat_norm1_records = vector('list',nj)
fhat_norm2_records = vector('list',nj)
fhat_norm3_records = vector('list',nj)
fhat_mm_records = vector('list',nj)
finf_records = matrix(0,nrow=nj,ncol=4)
f2_records = matrix(0,nrow=nj,ncol=4)
sig_records = matrix(0,nrow=nj,ncol=4)
RunTime_records = matrix(0,nrow=nj,ncol=4)

norm1_obj_records = vector('list',nj)
norm2_obj_records = vector('list',nj)
norm3_obj_records = vector('list',nj)
mm_obj_records = vector('list',nj)

## Experiments
for (j in j_lst){
  ## mcmc of norming dags
  Xj = matrix(X_norm[1:(2^j+1)],nrow=2^j+1,ncol=1)
  Yj = Y_norm[1:(2^j+1)]
  log_tau_fun <- function(tau,n,nu){
    return( -0.5 * n^(1/(2*nu+1))*tau^(2*nu*1/(2*nu+1)) )
  }

  dagj1 = dag_norm1_cpp[1:(2^j+1)]
  mcmc_obj1 = mcmc(Xj, dagj1, Yj, log_tau_fun, tau_bound_const, sig2_prior, sig2_bound, nu=nu1, tau=tau, sig=sig,
                   n_mcmc=n_mcmc, n_burn=n_burn)
  dagj2 = dag_norm2_cpp[1:(2^j+1)]
  mcmc_obj2 = mcmc(Xj, dagj2, Yj, log_tau_fun, tau_bound_const, sig2_prior, sig2_bound, nu=nu2, tau=tau, sig=sig,
                   n_mcmc=n_mcmc, n_burn=n_burn)
  dagj3 = dag_norm3_cpp[1:(2^j+1)]
  mcmc_obj3 = mcmc(Xj, dagj3, Yj, log_tau_fun, tau_bound_adapt, sig2_prior, sig2_bound, nu=nu3, tau=tau, sig=sig,
                   n_mcmc=n_mcmc, n_burn=n_burn)

  ## build maximin dag
  ordj_mm = GPvecchia::order_maxmin_exact(Xj)   ## ordj_mm converts order norming to order mm
  ordj_mm_2_ord_norming = sort(ordj_mm,index.return=TRUE)$ix
  mj = round(log(2^j+1)*2)
  Xj_mm = matrix(Xj[ordj_mm,],nrow=2^j+1,ncol=1)
  Yj_mm = Yj[ordj_mm]
  dagj_mm_mat = GPvecchia:::findOrderedNN_kdtree2(Xj_mm, mj)
  dagj_mm = vector('list',2^j+1)
  dagj_mm[[1]] = vector('numeric',0)
  for (i in 2:mj){
    dagj_mm[[i]] = dagj_mm_mat[i,2:i]
  }
  for (i in (mj+1):(2^j+1)){
    dagj_mm[[i]] = dagj_mm_mat[i,2:(mj+1)]
  }
  dagj_mm_cpp = Rdag_to_Cppdag(dagj_mm)

  mcmc_obj_mm = mcmc(Xj_mm, dagj_mm_cpp, Yj_mm, log_tau_fun, tau_bound_const, sig2_prior, sig2_bound, nu=beta, tau=tau, sig=sig, cov_type="matern",
                     n_mcmc=n_mcmc, n_burn=n_burn)

  ## records mcmc outputs
  j_ind = j-j_lst[1]+1
  fhat_norm1_records[[j_ind]] = confidence_bands(mcmc_obj1$Z_mcmc)$mean
  fhat_norm2_records[[j_ind]] = confidence_bands(mcmc_obj2$Z_mcmc)$mean
  fhat_norm3_records[[j_ind]] = confidence_bands(mcmc_obj3$Z_mcmc)$mean
  fhat_mm_records[[j_ind]] = confidence_bands(mcmc_obj_mm$Z_mcmc[,ordj_mm_2_ord_norming])$mean  ## fhat_mm is converted to ordering of Xj in this line

  sig_records[j_ind,] = c(mean(mcmc_obj1$sig_mcmc), mean(mcmc_obj2$sig_mcmc), mean(mcmc_obj3$sig_mcmc), mean(mcmc_obj_mm$sig_mcmc))
  finf_records[j_ind,] = c(max(abs(fhat_norm1_records[[j_ind]] - Yf_norm[1:(2^j+1)])),
                           max(abs(fhat_norm2_records[[j_ind]] - Yf_norm[1:(2^j+1)])),
                           max(abs(fhat_norm3_records[[j_ind]] - Yf_norm[1:(2^j+1)])),
                           max(abs(fhat_mm_records[[j_ind]] - Yf_norm[1:(2^j+1)])) )
  f2_records[j_ind,] = c(sqrt(sum((fhat_norm1_records[[j_ind]] - Yf_norm[1:(2^j+1)])^2)/(2^j+1)),
                         sqrt(sum((fhat_norm2_records[[j_ind]] - Yf_norm[1:(2^j+1)])^2)/(2^j+1)),
                         sqrt(sum((fhat_norm3_records[[j_ind]] - Yf_norm[1:(2^j+1)])^2)/(2^j+1)),
                         sqrt(sum((fhat_mm_records[[j_ind]] - Yf_norm[1:(2^j+1)])^2)/(2^j+1)) )
  RunTime_records[j_ind,] = c(mcmc_obj1$RunTime, mcmc_obj2$RunTime, mcmc_obj3$RunTime, mcmc_obj_mm$RunTime)

  norm1_obj_records[[j_ind]] = mcmc_obj1
  norm2_obj_records[[j_ind]] = mcmc_obj2
  norm3_obj_records[[j_ind]] = mcmc_obj3
  mm_obj_records[[j_ind]] = mcmc_obj_mm

  ## output for current j
  print(paste('n=',2^j+1))
  print(paste('l^2 estimation error:', f2_records[j_ind,]))
  print(paste('sigma, truth:',sig,' Estimates:', sig_records[j_ind,]))
  print(paste('Running Time:', RunTime_records[j_ind,]))
}

## plot estimation error and computation time
df_err = data.frame(logn=log(2^j_lst+1), BaseNorming=f2_records[,1], Oversmooth=f2_records[,2], Adaptation=f2_records[,3], BaseMaximin=f2_records[,4])
df_err = df_err[2:length(j_lst),]
df_err_melt = reshape2::melt(df_err, id.vars='logn', variable.name="Method", value.name="Error")
p_est =
  ggplot(df_err_melt) +
  geom_line(aes(x=logn,y=Error,color=Method)) +
  geom_point(aes(x=logn,y=Error,color=Method),shape=2) +
  theme_minimal(base_size = 20) +
  xlab('ln(n)') + ylab('Posterior Estimation Error') +
  theme(legend.position='none')

df_time = data.frame(logn=log(2^j_lst+1), BaseNorming=RunTime_records[,1], Oversmooth=RunTime_records[,2], Adaptation=RunTime_records[,3], BaseMaximin=RunTime_records[,4])
df_time = df_time[2:length(j_lst),]
df_time_melt = reshape2::melt(df_time, id.vars='logn', variable.name="Method", value.name="Time")
p_time =
  ggplot(df_time_melt) +
  geom_line(aes(x=logn,y=Time,color=Method)) +
  geom_point(aes(x=logn,y=Time,color=Method),shape=2) +
  theme_minimal(base_size = 20) +
  xlab('ln(n)') + ylab('Running Time')

grid.arrange(p_est, p_time, nrow=1, widths=c(1,1.32))

