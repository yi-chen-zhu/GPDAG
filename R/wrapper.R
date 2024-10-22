
#' DAG GP on grid
#'
#' This is the high-level function to perform DAG GP on grid, including both building the DAG and the MCMC on the DAG.
#'
#' @param Y numeric vector or numeric matrix, the response, with each element representing the value of response on its corresponding grid location.
#' @param minsep numeric value, the minimal separation distance of the grid.
#' @param adapt bool, whether adaptation on rescaling parameter is required. Default is FALSE.
#' @param nu numeric value, the smoothness of the Matern process
#' @param tau numeric value, initial value for rescaling parameter of the Matern process
#' @param sig numeric value, initial value for the standard deviation of nugget effects
#' @param tau_bound numeric vector of length 2, lower and upper bounds for tau
#' @param nug_prior numeric vector of length 2, parameters for inverse gamma prior for sigma^2
#' @param sig_bound numeric vector of length 2, lower and upper bounds for sigma, currently not used in the codes
#' @param n_mcmc unsigned integer, number of MCMC iterations
#' @param n_burn unsigned integer, number of MCMC iterations that are burnt
#' @param tol numeric value, the tolerance of squared l2 norm of conjugate gradient
#' @param ncgiter integer, the maximal number of conjugate gradient iterations
#' @param log_tau_prior_fun Function, log prior function for tau. Does not need to be specified if no adaptation is required.
#' @return Rcpp list consists results of MCMC:
#'         tau_mcmc, the MCMC samples of tau;
#'         sig_mcmc, the MCMC samples of sig;
#'         sig2_mean_mcmc, the mean value of sigma squared in each MCMC iteration;
#'         Z_mcmc, a matrix, each row contains the MCMC sample of the latent function in that iteration;
#'         Accept_ratio, average acceptance ratio for sampling the rescaling paramter tau. Is irrelevant if not in a hierarchical Bayes framework;
#'         L & D, the ouput of function DAG_Chol in the last MCMC iterations;
#'         RunTime: the total real world running time of the MCMC;
#'         CG_iters: the number of conjuate gradient iterations for each MCMC step;
#'         CG_err: the conjugate gradient error for each MCMC step;
#'         X: coordinates of the grid, ordered by their coordinate values;
#'        dag_ord: a list of numeric vectors, each entry contains the indices of its parent sets in DAG ordering;
#'        Y: the input response Y;
#'        sort_in_ord: a numeric vector, the ordering of DAG with respect to the coordinate ordering.
#' @useDynLib GPDAG, .registration = TRUE
#' @export
GPgrid <- function(Y, minsep, adapt=FALSE, nu=3/2, sig=0.1,
                    tau=10, tau_bound=NULL,
                    nug_prior=NULL, sig_bound=NULL,
                    n_mcmc=2000, n_burn=1000,
                    tol=1e-12, ncgiter=200,
                    log_tau_fun=NULL){
  if (!adapt){
    tau_bound = c(-1e-6,1e-6) + tau
    log_tau_fun <- function(tau){
      return(1)
    }
  }
  if (is.null(nug_prior)){
    sig_prior_size = 10
    nug_prior = c(sig_prior_size, sig_prior_size*sig^2)
  }
  if (is.null(sig_bound)){
    sig_bound = c(1e-6, sd(Y))
  }

  dag_obj = DAGgrid(Y, minsep, nu)
  X_ord = dag_obj$X_ord
  dag_ord = dag_obj$dag_ord
  Y_ord = dag_obj$Y_ord

  output = mcmc(matrix(X_ord), Rdag_to_Cppdag(dag_ord), Y_ord, log_tau_fun, tau_bound, nug_prior, sig_bound, nu=nu, tau=tau, sig=sig,
                n_mcmc=n_mcmc, n_burn=n_burn)
  output$X_ord = X_ord
  output$dag_ord = dag_ord
  output$Y_ord = Y_ord
  output$sort_in_ord = dag_obj$sort_in_ord

  return(output)
}


GPgeneral <- function(X, Y, adapt=FALSE, nu=3/2, sig=0.1,
                   tau=10, tau_bound=NULL,
                   nug_prior=NULL, sig_bound=NULL,
                   n_mcmc=4000, n_burn=2000,
                   tol=1e-12, ncgiter=200){
  if (!adapt){
    tau_bound = c(-1e-6,1e-6) + tau
  }
  if (is.null(nug_prior)){
    sig_prior_size = 10
    nug_prior = c(sig_prior_size, sig_prior_size*sig^2)
  }
  if (is.null(sig_bound)){
    sig_bound = c(1e-6, sd(Y))
  }

  dag_obj = DAGgeneral(X, nu)
  X_ord = dag_obj$X_ord
  dag_ord = dag_obj$dag_ord
  id2ord = dag_obj$ord
  Y_ord = Y[id2ord]
  ord2id = sort(id2ord,index.return=TRUE)$ix

  output = mcmc(X_ord, dag_ord, Y_ord, log_tau_fun, tau_bound, sig2_prior, sig2_bound, nu=nu, tau=tau, sig=sig,
                n_mcmc=n_mcmc, n_burn=n_burn)
  output$sort_in_ord = id2ord
  output$dag_ord = dag_ord
  output$X = X
  output$Y = Y
  output$Z_mcmc = output$Z_mcmc[,ord2id]

  return(output)
}


