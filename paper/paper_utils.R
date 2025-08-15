
##-----------------------utility functions for experiments---------------------------
confidence_bands <- function(Y, alpha=0.05){
  n_mcmc = nrow(Y)
  n = ncol(Y)
  i_low = max(floor(n_mcmc*alpha/2),1)
  i_up = floor(n_mcmc*(1-alpha/2))+1
  Ymean = rep(0,n)
  Ylow = rep(0,n)
  Yup = rep(0,n)
  for (i in 1:n){
    Y_now = sort(Y[,i])
    Ymean[i] = mean(Y_now)
    Ylow[i] = Y_now[i_low]
    Yup[i] = Y_now[i_up]
  }
  return(list(mean=Ymean,low=Ylow,up=Yup))
}

matrix_half <- function(M){
  eigenobj = eigen(M)
  U = eigenobj$vectors
  Dvec = pmax(eigenobj$values,0)
  Dhalf = diag(sqrt(Dvec))
  return( U %*% Dhalf %*% t(U))
}

matrix_tr<- function(M){
  return(sum(diag(M)))
}

matrix_W22 <- function(cov1, cov2){
  cov1_half = matrix_half(cov1)
  cov_mix = matrix_half( cov1_half %*% cov2 %*% cov1_half)
  return(matrix_tr(cov1) + matrix_tr(cov2) - 2*matrix_tr(cov_mix))
}
