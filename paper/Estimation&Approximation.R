library(Matrix)
library(ggplot2)
library(scales)
library(gridExtra)
library(GPDAG)
source("/home/yichen/GPDAG/paper/paper_utils.R")


## model and priors
nu = 1.5
beta = 1.5
sig = 0.1
sig_prior_size = 10
sig2_prior = c(sig_prior_size, sig_prior_size*sig^2)   ## gamma distribution parameters for sig2^{-1}
tau = 10
tau_bound = c(tau-1e-3,tau+1e-3)
log_tau_fun <- function(tau, n, nu){
  return(1)
}

## dataset; build norming dag
set.seed(1110)
J = 11
n = 2^J+1
dag_norming_obj = DAGgrid_per_1d(J, nu)
X = matrix(dag_norming_obj$X_ord,nrow=n,ncol=1)
ord_2_id = sort(X,index.return=TRUE)$ix
dag_norming = dag_norming_obj$dag_ord
dag_norming_cpp = Rdag_to_Cppdag(dag_norming)
X_cov = cov_matern(X, X, nu, tau)
Yf = MASS::mvrnorm(1,rep(0,n),X_cov)
Yf = Yf - mean(Yf)
Y = Yf + rnorm(n)*sig
sig2_bound = c(1e-3, sqrt(sum((Y-mean(Y))^2)/(n-1)) )

## Experimental settings
n_mcmc = 2000
n_burn = 1000
j_lst = seq(4,J)
nj = length(j_lst)
fhat_records = vector('list',nj)
finf_records = rep(0,nj)
f2_records = rep(0,nj)
sig_records = rep(0,nj)
fhat_mm_records = vector('list',nj)
finf_mm_records = rep(0,nj)
f2_mm_records = rep(0,nj)
sig_mm_records = rep(0,nj)
df_records = vector('list',nj)
plots = vector('list',nj)
RunTime_records = matrix(0,nrow=nj,ncol=2)
for (j in j_lst){
  ## mcmc of norming dag
  Xj = matrix(X[1:(2^j+1)],nrow=2^j+1,ncol=1)
  dagj = dag_norming_cpp[1:(2^j+1)]
  Yj = Y[1:(2^j+1)]
  mcmc_obj = mcmc(Xj, dagj, Yj, log_tau_fun, tau_bound, sig2_prior, sig2_bound, nu=nu, tau=tau, sig=sig,
                  n_mcmc=n_mcmc, n_burn=n_burn)
  fhat_obj = confidence_bands(mcmc_obj$Z_mcmc)
  sighat = mean(mcmc_obj$sig_mcmc)

  ## build maximin dag
  ordj_mm = GPvecchia::order_maxmin_exact(Xj)   ## ordjmm converts order norming to order mm
  ordj_mm_2_ord_norming = sort(ordj_mm,index.return=TRUE)$ix
  mj = round(log(2^j+1)*2)
  Xj_mm = matrix(Xj[ordj_mm,],nrow=2^j+1,ncol=1)
  Yj_mm = Yj[ordj_mm]
  # ord_mm_2_id = sort(X_mm,index.return=TRUE)$ix
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

  ## mcmc of maxmin dag
  mcmc_obj_mm = mcmc(Xj_mm, dagj_mm_cpp, Yj_mm, log_tau_fun, tau_bound, sig2_prior, sig2_bound, nu=nu, tau=tau, sig=sig,
                  n_mcmc=n_mcmc, n_burn=n_burn)
  fhat_mm_obj = confidence_bands(mcmc_obj_mm$Z_mcmc[,ordj_mm_2_ord_norming]) ## fhat_mm_obj is converted to ordering of Y in this line
  sighat_mm = mean(mcmc_obj_mm$sig_mcmc)

  ## records mcmc outputs
  j_ind = j-j_lst[1]+1
  fhat_records[[j_ind]] = fhat_obj$mean
  finf_records[j_ind] = max(abs(fhat_obj$mean-Yf[1:(2^j+1)]))
  f2_records[j_ind] = sqrt( sum((fhat_obj$mean-Yf[1:(2^j+1)])^2)/(2^j+1) )
  sig_records[j_ind] = sighat
  fhat_mm_records[[j_ind]] = fhat_mm_obj$mean
  finf_mm_records[j_ind] = max(abs(fhat_mm_obj$mean-Yf[1:(2^j+1)]))
  f2_mm_records[j_ind] = sqrt( sum((fhat_mm_obj$mean-Yf[1:(2^j+1)])^2)/(2^j+1) )
  sig_mm_records[j_ind] = sighat_mm
  RunTime_records[j_ind,] = c(mcmc_obj$RunTime, mcmc_obj_mm$RunTime)
  df <- data.frame(X=Xj, Y=Yj, Truth=Yf[1:(2^j+1)], Norming=fhat_obj$mean, Norming_low=fhat_obj$low, Norming_up=fhat_obj$up,
                   Maximin=fhat_mm_obj$mean, Maximin_low=fhat_mm_obj$low, Maximin_up=fhat_mm_obj$up)
  df_records[[j_ind]] = df

  ## output for current j
  print(paste('n=',2^j+1))
  print(paste('l^2 estimation error, Norming DAG:', f2_records[j_ind],
              ', Maximin DAG', f2_mm_records[j_ind] ))
  print(paste('sigma, truth:',sig,' Norming DAG:',sighat, ', Maximin DAG:',sighat_mm ))
  print(paste('Running Time, Norming DAG:',mcmc_obj$RunTime, ', Maximin DAG:',mcmc_obj_mm$RunTime ))
}

ylims = c(min(df_records[[2]][,2:ncol(df_records[[2]])])-0.1, max(df_records[[2]][,2:ncol(df_records[[1]])])+0.1)
# ylims = c(-1.3,1.3)
df_f = data.frame(X=X, Y=Yf)
for (j_ind in 1:length(j_lst)){
  df = df_records[[j_ind]]
  df_fits = reshape2::melt(df[c('X', 'Norming', 'Maximin')], id.vars = "X", variable.name = "Method", value.name = "Y")
  plots[[j_ind]] =
    ggplot() +
    geom_ribbon(data=df, mapping=aes(x=X,ymin=Norming_low,ymax=Norming_up),fill='#619CFF',alpha=0.4) +
    geom_ribbon(data=df, mapping=aes(x=X,ymin=Maximin_low,ymax=Maximin_up),fill='#F8766D',alpha=0.4) +
    geom_line(data=df_fits, mapping=aes(x=X,y=Y,color=Method), size=1) +
    geom_line(data=df_f, mapping=aes(x=X,y=Yf), color='black', size=1) +
    # geom_point(data=df, mapping=aes(x=X,y=Y), shape=20, size=1.2, alpha=1) +
    scale_colour_manual(values = c("#619CFF","#F8766D","black")) + ylim(ylims[1],ylims[2])+
    theme_minimal(base_size = 20) + labs(caption=paste('n=',2^j_lst[[j_ind]]+1)) +
    theme(legend.position='none', axis.title.x=element_blank(),axis.title.y=element_blank(),
          plot.caption=element_text(hjust=0.5, size=rel(1)))
}

df = df_records[[6]]
df_fits = reshape2::melt(df[c('X', 'Norming', 'Maximin', 'Truth')], id.vars = "X", variable.name = "Method", value.name = "Y")
p_jmax =
  ggplot() +
  geom_ribbon(data=df, mapping=aes(x=X,ymin=Norming_low,ymax=Norming_up),fill='#619CFF',alpha=0.4) +
  geom_ribbon(data=df, mapping=aes(x=X,ymin=Maximin_low,ymax=Maximin_up),fill='#F8766D',alpha=0.4) +
  geom_line(data=df_fits, mapping=aes(x=X,y=Y,color=Method), size=0.8) +
  scale_colour_manual(values = c("#619CFF","#F8766D","black")) + ylim(ylims[1],ylims[2])+
  theme_minimal(base_size = 20) + labs(caption=paste('n=',2^j_lst[[6]]+1)) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.caption=element_text(hjust=0.5, size=rel(1)))

grid.arrange(plots[[2]], plots[[4]], p_jmax, nrow=1, widths=c(1,1,1.23))

## plot estimation error and computation time
df_err = data.frame(logn=log(2^j_lst+1), Norming=f2_records, Maximin=f2_mm_records)
df_err = df_err[2:length(j_lst),]
df_err_melt = reshape2::melt(df_err, id.vars='logn', variable.name="Method", value.name="Error")
p_est =
  ggplot(df_err_melt) +
  geom_line(aes(x=logn,y=Error,color=Method)) +
  geom_point(aes(x=logn,y=Error,color=Method),shape=2) +
  theme_minimal(base_size = 20) + ylim(0,0.1) +
  xlab('ln(n)') + ylab('Posterior Estimation Error') +
  theme(legend.position='none')

df_time = data.frame(logn=log(2^j_lst+1), Norming=RunTime_records[,1], Maximin=RunTime_records[,2])
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
W22_table = matrix(0, nrow=length(jprior_lst), ncol=2)

for (j in jprior_lst){
  j_ind = j-jprior_lst[1]+1

  ## Norming dag
  Xj = X[1:(2^j+1)]
  dag_j = dag_norming[1:(2^j+1)]
  Xjtest = X[(2^j+2):n]
  if (j < J){
    dag_j_obj = DAGgrid_test_1d(Xj,Xjtest,nu)
    dag_j_cpp = Rdag_to_Cppdag(c(dag_j,dag_j_obj))
  } else{
    dag_j_cpp =  Rdag_to_Cppdag(dag_j)
  }
  chol_j_obj = DAG_Chol(X, dag_j_cpp, nu=nu, tau=tau)
  L_j = as.matrix(Matrix(chol_j_obj$L,sparse=FALSE))
  L_j_inv = solve(L_j)
  D_j = chol_j_obj$D
  cov_j = t(L_j_inv) %*% diag(D_j) %*% L_j_inv
  W22_table[j_ind,1] = matrix_W22(X_cov,cov_j)

  ## Maximin dag
  ordj_mm = GPvecchia::order_maxmin_exact(matrix(X[1:(2^j+1)],nrow=2^j+1,ncol=1))   ## ordjmm converts order norming to order mm
  ordj_mm_2_ord_norming = sort(ordj_mm,index.return=TRUE)$ix
  mj = round(log(2^j+1)*2)
  Xj_mm = matrix(X[ordj_mm],nrow=2^j+1,ncol=1)
  Yj_mm = Yj[ordj_mm]
  # ord_mm_2_id = sort(X_mm,index.return=TRUE)$ix
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
  W22_table[j_ind,2] = matrix_W22(X_cov_j_mm, cov_j_mm)
  print(paste('j=',j,' W22 of Norming:', W22_table[j_ind,1], ' W22 of Maximin:', W22_table[j_ind,2]))
}

W22_table_nonzero = pmax(W22_table, 0)
df_approx = data.frame(logn=log(2^jprior_lst+1), Norming=W22_table_nonzero[,1], Maximin=W22_table_nonzero[,2])
df_approx_melt = reshape2::melt(df_approx, id.vars='logn', variable.name="Method", value.name="W22")
p_approx =
  ggplot(df_approx_melt) +
  geom_line(aes(x=logn,y=W22,color=Method)) +
  geom_point(aes(x=logn,y=W22,color=Method),shape=2) +
  theme_minimal(base_size = 20) + ylim(0,30) +
  xlab('ln(n)') + ylab('Prior Approximation Error') +
  theme(legend.position='none')



grid.arrange(p_est, p_approx, p_time, nrow=1, widths=c(1,1,1.23))














