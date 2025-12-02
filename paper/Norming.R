library(Matrix)
library(ggplot2)
library(scales)
library(gridExtra)
library(latex2exp)
library(Rcpp)
library(GPDAG)
library(waveslim)
source("/home/yichen/GPDAG/paper/paper_utils.R")

## function to generate Holder-alpha functions
sample_wavelet_2d <- function(ncoord=129, J = 8, alpha = 1.5, wf = "d8") {
  if (ncoord> 2^J){
    J = floor(log2(ncoord)) + 1
  }
  N <- 2^J
  sup_adjust=N

  dummy <- matrix(0, nrow = N, ncol = N)
  dwt_struct <- dwt.2d(dummy, wf = wf, J = J)
  coeffs <- dwt_struct

  for (j in 1:J) {
    decay <- 2^(-(J-j) * (alpha + 1))  # controls Hölder α smoothness; already uses dimension d=2
    for (band in c("LH", "HL", "HH")) {
      nm <- paste0(band, j)
      sz <- dim(coeffs[[nm]])
      coeffs[[nm]] <- matrix(runif(prod(sz),min=-decay*sup_adjust, max=decay*sup_adjust), nrow = sz[1])
    }
  }
  coeffs[[paste0('LL',J)]] = matrix(runif(1, min=-sup_adjust, max=sup_adjust),nrow=1)

  f <- idwt.2d(coeffs)
  return(f[1:ncoord, 1:ncoord])
}

## model and priors
d = 2
nu = 2.5
beta = 2
l = ceiling(nu) - 1
m = choose(l+d,l)
sig = 0.1
sig_prior_size = 10
sig2_prior = c(sig_prior_size, sig_prior_size*sig^2)   ## gamma distribution parameters for sig2^{-1}
# sig2_bound = c(1e-3, sqrt(sum((Y-mean(Y))^2)/(n-1)) )
sig2_bound = c(1e-3, sig^2*100 )
tau = 8
tau_bound = c(tau-1e-3,tau+1e-3)
log_tau_fun <- function(tau, n, nu){
  return(1)
}

## dataset; build norming dag
set.seed(1110)
J = 8
n = (2^J+1)^d
dag_norming_obj = DAGgrid_per_hd(J, nu)
X = dag_norming_obj$X_ord
dag_norming = dag_norming_obj$dag_ord
dag_norming_cpp = Rdag_to_Cppdag(dag_norming)


## generate truth
Yf = rep(0, n)
for (i in 1:n){
  Yf[i] = X[i,2]-0.5 + 8*sin(X[i,1]*pi)
}
Y = Yf + rnorm(n)*sig

## Experimental settings
n_mcmc = 2000
n_burn = 1000
j_lst = seq(4,J)
nj = length(j_lst)
fhat_records = vector('list',nj)
finf_records = rep(0,nj)
f2_records = rep(0,nj)
sig_records = rep(0,nj)
fhat_nope_records = vector('list',nj)
finf_nope_records = rep(0,nj)
f2_nope_records = rep(0,nj)
sig_nope_records = rep(0,nj)
df_records = vector('list',nj)
plots = vector('list',nj)
tol=1e-16
maxcgiter=400
RunTime_records = matrix(0,nrow=nj,ncol=2)

for (j in j_lst){
  ## initialize dataset for jth experiment
  nj = (2^j+1)^d
  Xj = matrix(X[1:nj,],nrow=nj,ncol=d)
  Yj = Y[1:nj]
  Yfj = Yf[1:nj]

  ## mcmc of norming dag
  dagj = dag_norming_cpp[1:nj]
  mcmc_obj = mcmc(Xj, dagj, Yj, log_tau_fun, tau_bound, sig2_prior, sig2_bound, nu=nu, tau=tau, sig=sig,
                  n_mcmc=n_mcmc, n_burn=n_burn, tol=tol, maxcgiter=maxcgiter)
  fhat_obj = confidence_bands(mcmc_obj$Z_mcmc)
  sighat = mean(mcmc_obj$sig_mcmc)

  ## mcmc of none-norming dag
  dagj_nope = Rdag_to_Cppdag(DAGgrid_Not_Norm(j)$dag_ord)
  mcmc_nope_obj = mcmc(Xj, dagj_nope, Yj, log_tau_fun, tau_bound, sig2_prior, sig2_bound, nu=nu, tau=tau, sig=sig,
                  n_mcmc=n_mcmc, n_burn=n_burn, tol=tol, maxcgiter=maxcgiter)
  fhat_nope_obj = confidence_bands(mcmc_nope_obj$Z_mcmc)
  sighat_nope = mean(mcmc_nope_obj$sig_mcmc)

  ## records mcmc outputs
  j_ind = j-j_lst[1]+1
  fhat_records[[j_ind]] = fhat_obj$mean
  finf_records[j_ind] = max(abs(fhat_obj$mean-Yfj))
  f2_records[j_ind] = sqrt( sum((fhat_obj$mean-Yfj)^2)/length(Yfj) )
  sig_records[j_ind] = sighat
  fhat_nope_records[[j_ind]] = fhat_nope_obj$mean
  finf_nope_records[j_ind] = max(abs(fhat_nope_obj$mean-Yfj))
  f2_nope_records[j_ind] = sqrt( sum((fhat_nope_obj$mean-Yfj)^2)/length(Yfj) )
  sig_nope_records[j_ind] = sighat_nope
  RunTime_records[j_ind,] = c(mcmc_obj$RunTime, mcmc_nope_obj$RunTime)
  df <- data.frame(X=Xj, Y=Yj, Truth=Yfj, Norming=fhat_obj$mean, Norming_low=fhat_obj$low, Norming_up=fhat_obj$up,
                   Maximin=fhat_nope_obj$mean, Maximin_low=fhat_nope_obj$low, Maximin_up=fhat_nope_obj$up)
  df_records[[j_ind]] = df

  ## output for current j
  print(paste('sqrt(n)=',2^j+1))
  print(paste('l^2 estimation error, Norming DAG:', f2_records[j_ind],
              ', N-N DAG', f2_nope_records[j_ind] ))
  print(paste('sigma, truth:',sig,' Norming DAG:',sighat, ', N-N DAG:',sighat_nope ))
  print(paste('Running Time, Norming DAG:',mcmc_obj$RunTime, ', N-N DAG:',mcmc_nope_obj$RunTime ))
}

plot_size = 1.3
df1 = data.frame(logn=2*log(2^j_lst+1), Norming=f2_records, Rect=f2_nope_records)
df1 = df1[2:nrow(df1),]
df1_melt = reshape2::melt(df1, id.vars='logn', variable.name="Method", value.name="Error")
p1 =
  ggplot(df1_melt) +
  geom_line(aes(x=logn,y=Error,color=Method),size=1) +
  geom_point(aes(x=logn,y=Error,color=Method),shape=2,size=2.5) +
  theme_minimal(base_size = 20) +
  xlab('ln(n)') + ylab('Posterior Estimation Error') +
  ylim(0,max(df1_melt[,3])) +
  theme(legend.position='none')







##---------------------test minimal eigenvalues------------------------------
Xt1 = t(matrix(c(0,0,0,1,1,0,1,1,2,0,2,1),nrow=2))
Xt2 = t(matrix(c(0,0,0,1,1,0,1,1,2,1,1,2),nrow=2))
j_lst = c(5,6,7,8,9,10,11,12)
n_scales = length(j_lst)
eigen_mat = matrix(0, nrow=n_scales, ncol=2)
for (i in 1:n_scales){
  j = j_lst[i]
  tau_j = tau/2^j
  covt1 = cov_matern(Xt1,Xt1,nu=5/2,tau=tau_j)
  covt2 = cov_matern(Xt2,Xt2,nu=5/2,tau=tau_j)
  eigen_mat[i,] = c(min(eigen(covt1)$values),min(eigen(covt2)$values))
}
logeigen_mat = log(eigen_mat)
df2 = data.frame(logn=2*log(2^j_lst+1), Norming=logeigen_mat[,2], Rect=logeigen_mat[,1])
df2_melt = reshape2::melt(df2, id.vars='logn', variable.name="Method", value.name="Eigen")
p2 =
  ggplot(df2_melt) +
  geom_line(aes(x=logn,y=Eigen,color=Method),size=1) +
  geom_point(aes(x=logn,y=Eigen,color=Method),shape=2,size=2.5) +
  theme_minimal(base_size = 20) +
  xlab('ln(n)') + ylab('Log Minimal Eigenvalue')


## plot examples of norming and rectangle(non-norming) sets
X_norming = t(matrix(
  c(0,0, 0,0.5, 0.5,0, 0.5,0.5, 0.5,1, 1,0.5), nrow=2
))
X_rect = t(matrix(
  c(0,0, 0,0.5, 0.5,0, 0.5,0.5, 1,0, 1,0.5), nrow=2
))


df = data.frame(x1=X_norming[,1],x2=X_norming[,2])
p_norming =
  ggplot(df) +
  geom_point(aes(x=x1,y=x2),color='red',size=4) +
  xlim(-0.02, 1.02) + ylim(-0.02, 1.02) +
  labs(caption=TeX(r'(Norming set: $c_N=7$)')) +
  theme_minimal(base_size = 20) +
  theme(legend.position='none',axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.caption=element_text(hjust=0.5, size=rel(1)))

df = data.frame(x1=X_rect[,1],x2=X_rect[,2])
p_rect =
  ggplot(df) +
  geom_point(aes(x=x1,y=x2),color='red',size=4) +
  xlim(-0.02, 1.02) + ylim(-0.02, 1.02) +
  labs(caption=TeX(r'(Rectangular set: $c_N=\infty$)')) +
  theme_minimal(base_size = 20) +
  theme(legend.position='none',axis.title.x=element_blank(),axis.title.y=element_blank(),
        plot.caption=element_text(hjust=0.5, size=rel(1)))

p_sets = grid.arrange(p_norming, p_rect, nrow=2)



grid.arrange(grid::nullGrob(), p_sets, grid::nullGrob(), p1, grid::nullGrob(), p2, nrow=1,
             widths=c(0.05, 0.5, 0.1, 1, 0.05, 1.23))










