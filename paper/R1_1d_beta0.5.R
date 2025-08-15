library(Matrix)
library(ggplot2)
library(scales)
library(gridExtra)
library(GPDAG)
source("/home/yichen/GPDAG/paper/paper_utils.R")

##------------------------ truth is Sobolev 1/2 smooth -------------------------
## generate the underlying truth
beta = 0.5
set.seed(77)
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
tau_bound = c(tau-1e-3,tau+1e-3)
log_tau_fun <- function(tau, n, nu){
  return(1)
}
sig2_bound = c(1e-3, sqrt(sum((Y-mean(Y))^2)/(n-1)) )

## Norming methods. Notice two norming methods share the same layer structures, but with different DAGs
nu1 = 0.5  ## result parent set size 1
nu2 = 1.5  ## result parent set size 2
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

## Experimental settings
n_mcmc = 2000
n_burn = 1000
j_lst = seq(4,J)
nj = length(j_lst)

fhat_norm1_records = vector('list',nj)
fhat_norm2_records = vector('list',nj)
fhat_mm_records = vector('list',nj)
finf_records = matrix(0,nrow=nj,ncol=3)
f2_records = matrix(0,nrow=nj,ncol=3)
sig_records = matrix(0,nrow=nj,ncol=3)
RunTime_records = matrix(0,nrow=nj,ncol=3)

norm1_obj_records = vector('list',nj)
norm2_obj_records = vector('list',nj)
mm_obj_records = vector('list',nj)

## Experiments
for (j in j_lst){
  ## mcmc of norming dags
  Xj = matrix(X_norm[1:(2^j+1)],nrow=2^j+1,ncol=1)
  Yj = Y_norm[1:(2^j+1)]

  dagj1 = dag_norm1_cpp[1:(2^j+1)]
  mcmc_obj1 = mcmc(Xj, dagj1, Yj, log_tau_fun, tau_bound, sig2_prior, sig2_bound, nu=beta, tau=tau, sig=sig,
                  n_mcmc=n_mcmc, n_burn=n_burn)
  dagj2 = dag_norm2_cpp[1:(2^j+1)]
  mcmc_obj2 = mcmc(Xj, dagj2, Yj, log_tau_fun, tau_bound, sig2_prior, sig2_bound, nu=beta, tau=tau, sig=sig,
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

  mcmc_obj_mm = mcmc(Xj_mm, dagj_mm_cpp, Yj_mm, log_tau_fun, tau_bound, sig2_prior, sig2_bound, nu=beta, tau=tau, sig=sig,
                  n_mcmc=n_mcmc, n_burn=n_burn)

  ## records mcmc outputs
  j_ind = j-j_lst[1]+1
  fhat_norm1_records[[j_ind]] = confidence_bands(mcmc_obj1$Z_mcmc)$mean
  fhat_norm2_records[[j_ind]] = confidence_bands(mcmc_obj2$Z_mcmc)$mean
  fhat_mm_records[[j_ind]] = confidence_bands(mcmc_obj_mm$Z_mcmc[,ordj_mm_2_ord_norming])$mean  ## fhat_mm is converted to ordering of Xj in this line

  sig_records[j_ind,] = c(mean(mcmc_obj1$sig_mcmc), mean(mcmc_obj2$sig_mcmc), mean(mcmc_obj_mm$sig_mcmc))
  finf_records[j_ind,] = c(max(abs(fhat_norm1_records[[j_ind]] - Yf_norm[1:(2^j+1)])),
                           max(abs(fhat_norm2_records[[j_ind]] - Yf_norm[1:(2^j+1)])),
                           max(abs(fhat_mm_records[[j_ind]] - Yf_norm[1:(2^j+1)])) )
  f2_records[j_ind,] = c(sqrt(sum((fhat_norm1_records[[j_ind]] - Yf_norm[1:(2^j+1)])^2)/(2^j+1)),
                         sqrt(sum((fhat_norm2_records[[j_ind]] - Yf_norm[1:(2^j+1)])^2)/(2^j+1)),
                         sqrt(sum((fhat_mm_records[[j_ind]] - Yf_norm[1:(2^j+1)])^2)/(2^j+1)) )
  RunTime_records[j_ind,] = c(mcmc_obj1$RunTime, mcmc_obj2$RunTime, mcmc_obj_mm$RunTime)

  norm1_obj_records[[j_ind]] = mcmc_obj1
  norm2_obj_records[[j_ind]] = mcmc_obj2
  mm_obj_records[[j_ind]] = mcmc_obj_mm

  ## output for current j
  print(paste('n=',2^j+1))
  print(paste('l^2 estimation error:', f2_records[j_ind,]))
  print(paste('sigma, truth:',sig,' Estimates:', sig_records[j_ind,]))
  print(paste('Running Time:', RunTime_records[j_ind,]))
}

## plot estimation error and computation time
df_err = data.frame(logn=log(2^j_lst+1), Norming1=f2_records[,1], Norming2=f2_records[,2], Maximin=f2_records[,3])
df_err = df_err[2:length(j_lst),]
df_err_melt = reshape2::melt(df_err, id.vars='logn', variable.name="Method", value.name="Error")
p_est =
  ggplot(df_err_melt) +
  geom_line(aes(x=logn,y=Error,color=Method)) +
  geom_point(aes(x=logn,y=Error,color=Method),shape=2) +
  theme_minimal(base_size = 20) +
  ylim(0,0.12) +
  xlab('ln(n)') + ylab('Posterior Estimation Error') +
  theme(legend.position='none')

df_time = data.frame(logn=log(2^j_lst+1), Norming1=RunTime_records[,1], Norming2=RunTime_records[,2], Maximin=RunTime_records[,3])
df_time = df_time[2:length(j_lst),]
df_time_melt = reshape2::melt(df_time, id.vars='logn', variable.name="Method", value.name="Time")
p_time =
  ggplot(df_time_melt) +
  geom_line(aes(x=logn,y=Time,color=Method)) +
  geom_point(aes(x=logn,y=Time,color=Method),shape=2) +
  theme_minimal(base_size = 20) +
  xlab('ln(n)') + ylab('Running Time')

##----------------------------prior approximation-------------------------------
jprior_lst = seq(5,J)
# X_cov_id = X_cov[ord_2_id,ord_2_id]
W22_table = matrix(0, nrow=length(jprior_lst), ncol=3)

for (j in jprior_lst){
  j_ind = j-jprior_lst[1]+1

  ## Norming dag 1
  Xj = X_norm[1:(2^j+1)]
  dag_j = dag_norm1[1:(2^j+1)]
  Xjtest = X_norm[(2^j+2):n]
  if (j < J){
    dag_j_obj = DAGgrid_test_1d(Xj,Xjtest,beta)
    dag_j_cpp = Rdag_to_Cppdag(c(dag_j,dag_j_obj))
  } else{
    dag_j_cpp =  Rdag_to_Cppdag(dag_j)
  }
  chol_j_obj = DAG_Chol(X_norm, dag_j_cpp, nu=beta, tau=tau)
  L_j = as.matrix(Matrix(chol_j_obj$L,sparse=FALSE))
  L_j_inv = solve(L_j)
  D_j = chol_j_obj$D
  cov_j = t(L_j_inv) %*% diag(D_j) %*% L_j_inv
  W22_table[j_ind,1] = matrix_W22(X_cov,cov_j)

  ## Norming dag 2
  Xj = X_norm[1:(2^j+1)]
  dag_j = dag_norm2[1:(2^j+1)]
  Xjtest = X_norm[(2^j+2):n]
  if (j < J){
    dag_j_obj = DAGgrid_test_1d(Xj,Xjtest,beta)
    dag_j_cpp = Rdag_to_Cppdag(c(dag_j,dag_j_obj))
  } else{
    dag_j_cpp =  Rdag_to_Cppdag(dag_j)
  }
  chol_j_obj = DAG_Chol(X_norm, dag_j_cpp, nu=beta, tau=tau)
  L_j = as.matrix(Matrix(chol_j_obj$L,sparse=FALSE))
  L_j_inv = solve(L_j)
  D_j = chol_j_obj$D
  cov_j = t(L_j_inv) %*% diag(D_j) %*% L_j_inv
  W22_table[j_ind,2] = matrix_W22(X_cov,cov_j)

  ## Maximin dag
  ordj_mm = GPvecchia::order_maxmin_exact(matrix(X_norm[1:(2^j+1)],nrow=2^j+1,ncol=1))   ## ordjmm converts order norming to order mm
  ordj_mm_2_ord_norming = sort(ordj_mm,index.return=TRUE)$ix
  mj = round(log(2^j+1)*2)
  Xj_mm = matrix(X_norm[ordj_mm],nrow=2^j+1,ncol=1)
  dagj_mm_mat = GPvecchia:::findOrderedNN_kdtree2(Xj_mm, mj)
  dagj_mm = vector('list',n)
  dagj_mm[[1]] = vector('numeric',0)
  for (i in 2:mj){
    dagj_mm[[i]] = dagj_mm_mat[i,2:i]
  }
  for (i in (mj+1):(2^j+1)){
    dagj_mm[[i]] = dagj_mm_mat[i,2:(mj+1)]
  }
  if (j<J){
    nn_mat = FNN::get.knnx(Xj_mm, Xjtest, mj)$nn.index
    for (i in (2^j+2):n){
      dagj_mm[[i]] = nn_mat[i-(2^j+1),]
    }
    Xj_mm_all = as.matrix(c(Xj_mm, Xjtest))
  } else{
    Xj_mm_all = Xj_mm
  }
  dagj_mm_cpp = Rdag_to_Cppdag(dagj_mm)
  chol_j_mm_obj = DAG_Chol(Xj_mm_all, dagj_mm_cpp, nu=nu, tau=tau)
  L_j_mm = as.matrix(Matrix(chol_j_mm_obj$L,sparse=FALSE))
  L_j_mm_inv = solve(L_j_mm)
  D_j_mm = chol_j_mm_obj$D
  cov_j_mm = t(L_j_mm_inv) %*% diag(D_j_mm) %*% L_j_mm_inv
  X_cov_j_mm = cov_matern(Xj_mm_all,Xj_mm_all,nu=nu,tau=tau)
  W22_table[j_ind,3] = matrix_W22(X_cov_j_mm, cov_j_mm)
  print(paste('j=',j,' W22 of Norming 1:', W22_table[j_ind,1], ' W22 of Norming 2:', W22_table[j_ind,2], ' W22 of Maximin:', W22_table[j_ind,3]))
}

W22_table_nonzero = pmax(W22_table, 0)
df_approx = data.frame(logn=log(2^jprior_lst+1), Norming1=W22_table_nonzero[,1], Norming2=W22_table_nonzero[,2], Maximin=W22_table_nonzero[,3])
df_approx_melt = reshape2::melt(df_approx, id.vars='logn', variable.name="Method", value.name="W22")
p_approx =
  ggplot(df_approx_melt) +
  geom_line(aes(x=logn,y=W22,color=Method)) +
  geom_point(aes(x=logn,y=W22,color=Method),shape=2) +
  theme_minimal(base_size = 20) +
  xlab('ln(n)') + ylab('Prior Approximation Error') +
  theme(legend.position='none')

grid.arrange(p_est, p_approx, p_time, nrow=1, widths=c(1,1,1.23))




