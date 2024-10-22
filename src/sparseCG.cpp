#include "RcppArmadillo.h"
#include <RcppEigen.h>
#include <cmath>
#include <tuple>
#include <iostream>
#include <boost/math/special_functions/bessel.hpp>
#include <chrono>  // for time measurement
using namespace Rcpp;

inline Eigen::VectorXd vec_arma_2_Eigen(arma::vec v){
  Eigen::VectorXd Eigen_v = Eigen::Map<Eigen::VectorXd>(v.memptr(), v.n_elem);
  return Eigen_v;
}


inline double scalar_Rfun(SEXP x){
  NumericVector v(x);
  return(v[0]);
}

//' Covariance matrix of Matern GP
//'
//' This function evaluates the covariance matrix of Matern GP using the coordinates.
//'
//' @param X1 arma::mat, each row representing the coordiantes of a location
//' @param X2 arma::mat, each row representing the coordiantes of a location
//' @param nu double, the smoothness of the Matern process
//' @param tau double, the time rescaling parameter of the Matern process
//' @param s double, the space rescaling parameter of the Matern process. Currently set as 1
//' @return the computed covariance matrix.
//' @useDynLib GPDAG, .registration = TRUE
//' @export
// [[Rcpp::export]]
arma::mat cov_matern(const arma::mat& X1,
                     const arma::mat& X2,
                     double nu=3/2,
                     double tau=1,
                     double s=1){
  arma::mat res(X1.n_rows,X2.n_rows);
  for(unsigned int i=0; i<X1.n_rows; i++){
    arma::rowvec cri = X1.row(i);
    for(unsigned int j=0; j<X2.n_rows; j++){
      arma::rowvec delta = cri - X2.row(j);
      double delta_x = arma::norm(delta) * tau;
      if(delta_x > 0.0){
        res(i, j) = s * std::pow(delta_x, nu) * std::pow(2,1-nu) / std::tgamma(nu) * boost::math::cyl_bessel_k(nu, delta_x);
      } else {
        res(i, j) = s;
      }
    }
  }
  return res;
}


//' Cholesky of DAG GP covariance matrix
//'
//' This function performs Cholesky decomposition of DAG GP covariance function. No input of the whole covariance matrix is
//' needed. Instead, this function will call evaluate the covariance matrix based on the DAG structure.
//'
//' @param X arma::mat, with each row representing the coordiante of a location
//' @param dag arma::field<arma::uvec>, each entry contains the indices of its parent sets in DAG ordering
//' @param nu double, the smoothness of the Matern process
//' @param tau double, rescaling parameter of the Matern process
//' @return Rcpp list consists of two entries, L and D, such that L D^{-1} L^T is the precision matrix.
//' D is the conditonal variance for ith element, while the off-diagonal elements of L[,i] are negative conditional regression coefficients of ith elements on its parents.
//' @useDynLib GPDAG, .registration = TRUE
//' @export
// [[Rcpp::export]]
Rcpp::List DAG_Chol(const arma::mat& X,
                    const arma::field<arma::uvec>& dag,
                    double nu=3/2,
                    double tau=1){
  unsigned int n = X.n_rows;
  Eigen::VectorXd D(n);
  unsigned int n_nonzero=n;
  for(int i=0; i<n; i++){
    n_nonzero += dag(i).n_elem;
  }
  Eigen::SparseMatrix<double> L(n,n);
  std::vector<Eigen::Triplet<double>> triples_L(n_nonzero);
// #ifdef _OPENMP
//   #pragma omp parallel for
// #endif
  for(int i=0; i<n; i++){
    const arma::mat& xi = X.row(i);
    const arma::mat& xpa = X.rows(dag(i));

    arma::vec K_pai = cov_matern(xpa, xi, nu, tau);
    arma::mat K_papa = cov_matern(xpa, xpa, nu, tau);
    arma::mat Kinv_papa = arma::inv_sympd(K_papa);
    arma::vec interp_coef = Kinv_papa * K_pai;
    double K_cond = arma::as_scalar(cov_matern(xi,xi) - K_pai.t() * interp_coef);
    D(i) = K_cond;
    double sqrtDi = sqrt(D(i));
    triples_L.push_back(Eigen::Triplet<double>(i,i,1));
    for(int j=0;j<dag(i).n_elem; j++){
      triples_L.push_back(Eigen::Triplet<double>(dag(i)(j),i,-interp_coef(j)));
    }
  }

  L.setFromTriplets(triples_L.begin(), triples_L.end());

  return List::create(
    Rcpp::Named("L") = L,
    Rcpp::Named("D") = D
  );
}


// solve linear systems: (I+D^{1/2}L^{-1}L^{-T}D^{1/2}/sig^2) D^{-1/2}L^TZ = W, return Z
// treat D^{-1/2}L^TZ as a whole during CG iterations; solve Z after CG converges
std::tuple<Eigen::VectorXd, unsigned int, double> DAG_CG(const Eigen::SparseMatrix<double>& L,
                       const Eigen::VectorXd& D,
                       const double sig,
                       const Eigen::VectorXd& W,
                       const double tol=1e-12,
                       const double maxiter=200){
  unsigned int n = L.rows();
  const Eigen::SparseMatrix<double> LT = L.transpose();
  const Eigen::VectorXd Dsqrt = D.array().sqrt();
  const double sig2inv = 1/(sig*sig);
  Eigen::VectorXd r = W;
  Eigen::VectorXd p = W;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  double r22 = r.squaredNorm();
  const double r220 = r22;
  double alpha = 0;
  double beta = 0;
  double tol_new = 1;

  unsigned int niter = 0;
  while(niter<maxiter){
    Eigen::VectorXd Bp = LT.triangularView<Eigen::Lower>().solve( Dsqrt.cwiseProduct(p) );
    Eigen::VectorXd BBp = Dsqrt.cwiseProduct( L.triangularView<Eigen::Upper>().solve(Bp) );
    Eigen::VectorXd Ap = p + BBp * sig2inv;
    double pAp = p.dot(Ap);
    alpha = r22 / pAp;
    x = x + alpha * p;
    Eigen::VectorXd r_new = r - alpha * Ap;
    double r22_new = r_new.squaredNorm();
    tol_new = r22_new/r220;
    if (tol_new < tol){
      break;
    }
    beta = r22_new / r22;
    r = r_new;
    r22 = r22_new;
    p = r + beta * p;
    niter += 1;
  }
  Eigen::VectorXd Z = LT.triangularView<Eigen::Lower>().solve( Dsqrt.cwiseProduct(x) );
  return std::make_tuple(Z, niter, tol_new);
}


//' MCMC for DAG GP
//'
//' This is the advanced function to perform MCMC for DAG GP and works for generic training set and DAG structures. For implementation on grid, one can call GPgrid instead.
//' This functions takes the coordinates of training locations (X), the DAG (dag) and the response (Y) as input. It performs MCMC and return
//' the MCMC samples as output. It uses a Gibbs sampler framework, whiling calling preconditioner Conjugate gradient when sampling latent process values.
//'
//' @param X arma::mat, with each row representing the coordiante of a location
//' @param dag arma::field<arma::uvec>, each entry contains the indices of its parent sets in DAG ordering
//' @param log_tau_prior_fun Function,log prior function for tau
//' @param tau_bound arma::vec, lower and upper bounds for tau
//' @param nug_prior arma::vec, parameters for inverse gamma prior for sigma^2
//' @param sig_bound arma::vec, lower and upper bounds for sigma, currently not used in the codes
//' @param nu double, the smoothness of the Matern process
//' @param tau double, initial value for rescaling parameter of the Matern process
//' @param sig double, initial value for the standard deviation of nugget effects
//' @param n_mcmc unsigned integer, number of MCMC iterations
//' @param n_burn unsigned integer, number of MCMC iterations that are burnt
//' @param tol double, the tolerance of squared l2 norm of conjugate gradient
//' @param maxcgiter unsigned int, the maximal number of conjugate gradient iterations
//' @return Rcpp list consists results of MCMC:
//'         tau_mcmc, the MCMC samples of tau;
//'         sig_mcmc, the MCMC samples of sig;
//'         sig2_mean_mcmc, the mean value of sigma squared in each MCMC iteration;
//'         Z_mcmc, a matrix, each row contains the MCMC sample of the latent function in that iteration;
//'         Accept_ratio, average acceptance ratio for sampling the rescaling paramter tau. Is irrelevant if not in a hierarchical Bayes framework;
//'         L & D, the ouput of function DAG_Chol in the last MCMC iterations;
//'         RunTime: the total real world running time of the MCMC;
//'         CG_iters: the number of conjuate gradient iterations for each MCMC step;
//'         CG_err: the conjugate gradient error for each MCMC step.
//' @useDynLib GPDAG, .registration = TRUE
//' @export
// [[Rcpp::export]]
Rcpp::List mcmc(const arma::mat& X,
                const arma::field<arma::uvec>& dag,
                const arma::vec& Y,
                Function log_tau_prior_fun,  // log prior function for tau
                const arma::vec tau_bound, // lower and upper bounds for tau as adaptation requires lower bound tau
                const arma::vec nug_prior,  // parameters of inverse gamma prior for sigma^2
                const arma::vec sig_bound,
                double nu=3/2,
                double tau=1, // initial value for tau
                double sig=0.1,  // initial value for sigma
                unsigned int n_mcmc=4000,  // default number of mcmc iterations
                unsigned int n_burn=2000,  // default number of burn-in iterations for mcmc
                const double tol=1e-12,
                const unsigned int maxcgiter=200
                ){
  auto t0 = std::chrono::high_resolution_clock::now();
  // data initialization
  unsigned int n = Y.n_elem;
  Eigen::VectorXd Ye = vec_arma_2_Eigen(Y);
  std::vector<Eigen::Triplet<double>> triple_id;
  for (int i=0; i<n; i++){
    triple_id.emplace_back(i,i,1);
  }
  Eigen::SparseMatrix<double> Id(n,n);
  Id.setFromTriplets(triple_id.begin(),triple_id.end());
  Eigen::VectorXd Z = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd Z_diff = Ye-Z;
  std::tuple<Eigen::VectorXd, unsigned int, double> cg_obj;

  // tau related variables initialization
  Rcpp::List chol_obj = DAG_Chol(X, dag, nu, tau);
  Eigen::SparseMatrix<double> L = chol_obj["L"];
  Eigen::VectorXd D = chol_obj["D"];
  Eigen::SparseMatrix<double> Phi = L * L.transpose();
  Eigen::SparseMatrix<double> Phi_nug = Phi + Id/(sig*sig);
  double log_detL = 0;
  for (int i=0; i<n; i++){
    log_detL += log(D(i));
  }
  double log_tau_prior = scalar_Rfun(log_tau_prior_fun(tau));

  // MCMC initialization
  unsigned int n_keep = n_mcmc - n_burn;
  unsigned int n_accept = 0;
  arma::vec tau_records(n_keep);
  arma::vec sig_records(n_keep);
  arma::vec sig2_mean_records(n_keep);
  Eigen::MatrixXd Z_records(n_keep,n);
  double AM_theta=1;  // initial log variance of proposal distribution of sig
  double AM_step0=0.5; // initial step size for adapatation of theta
  double AM_step=AM_step;
  double AM_alpha_star = 0.44; // asymptotically optimal acceptance rate for d=1
  Eigen::VectorXd W(n);
  Eigen::VectorXd W_precond(n);
  double logden;
  arma::vec cg_iters_records(n_mcmc);
  arma::vec cg_err_records(n_mcmc);

  // MCMC loops
  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> deltat = t1-t0;
  std::cout<<"head time: "<< deltat.count()<<std::endl;
  for (int imc=0; imc<n_mcmc; imc++){
    // sample Z with CG
    Eigen::VectorXd W1 = vec_arma_2_Eigen(arma::randn(n));
    Eigen::VectorXd W2 = vec_arma_2_Eigen(arma::randn(n));
    if (D.minCoeff()>sig*sig){
      // directly perform CG without preconditioners
      W = Ye/(sig*sig) + L*W1 + W2/sig;
      Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> cg;
      cg.setTolerance(sqrt(tol)/(sig*sig));
      cg.compute(Phi_nug);
      Z = cg.solveWithGuess(W,Ye);
      cg_iters_records(imc) = cg.iterations();
      cg_err_records(imc) = cg.error();
    } else{
      // perform CG with a preconditioner as the Cholesky decomposition of DAG covariance
      Eigen::VectorXd Dsqrt = D.array().sqrt();
      Eigen::VectorXd W_partial = Dsqrt.cwiseProduct( L.triangularView<Eigen::Upper>().solve(Ye/(sig*sig) + W2/sig) );
      W_precond = W_partial + W1;
      cg_obj = DAG_CG(L,D,sig,W_precond,tol,maxcgiter);
      Z = std::get<0>(cg_obj);
      cg_iters_records(imc) = std::get<1>(cg_obj);
      cg_err_records(imc) = std::get<2>(cg_obj);
    }

    Z_diff = Ye - Z;
    Eigen::MatrixXd LTZ = L.transpose() * Z;
    double ZPhiZ = LTZ.squaredNorm();
    logden = log_tau_prior + log_detL - ZPhiZ/2;

    // sample sig
    double sig2_mean = ( nug_prior(1)+Z_diff.squaredNorm()/2 )/(nug_prior(0)+n/2-1);
    sig = sqrt(1/R::rgamma(nug_prior(0)+n/2, 1/(nug_prior(1)+Z_diff.squaredNorm()/2)) );

    // sample tau
    double tau_prop = tau + exp(AM_theta)*arma::randn(1)(0);
    AM_step = AM_step0 / pow(imc+1,3/4);
    bool accept=false;
    if (tau_prop>tau_bound(0) and tau_prop<tau_bound(1)){
      Rcpp::List chol_obj_prop = DAG_Chol(X, dag, nu, tau_prop);
      Eigen::SparseMatrix<double> L_prop = chol_obj_prop["L"];
      Eigen::VectorXd D_prop = chol_obj_prop["D"];
      double log_detL_prop = 0;
      for (int i=0; i<n; i++){
        log_detL_prop += log(D_prop(i));
      }
      Eigen::MatrixXd LTZ_prop = L_prop.transpose() * Z;
      double ZPhiZ_prop = LTZ_prop.squaredNorm();
      double log_tau_prior_prop = scalar_Rfun(log_tau_prior_fun(tau_prop));
      double logden_prop = log_tau_prior_prop + log_detL_prop - ZPhiZ_prop/2;
      double logden_diff = logden_prop - logden;
      double AM_alpha = AM_alpha_star;
      if (logden_diff>0){
        accept=true;
        AM_alpha = 1;
      } else {
        AM_alpha = exp(logden_diff);
        accept= arma::randu<arma::vec>(1)(0) < AM_alpha;
      }
      if (accept){
        tau = tau_prop;
        L = L_prop;
        D = D_prop;
        Phi = L * L.transpose();
        Phi_nug = Phi + Id/(sig*sig);
        log_tau_prior = log_tau_prior_prop;
        logden = logden_prop;
      }
      AM_theta = AM_theta + (AM_alpha-AM_alpha_star)*AM_step;
    }

    // record tau, sig, Z
    if (imc>=n_burn){
      if (accept){
        n_accept += 1;
      }
      tau_records(imc-n_burn) = tau;
      sig_records(imc-n_burn) = sig;
      sig2_mean_records(imc-n_burn) = sig2_mean;
      Z_records.row(imc-n_burn) = Z.transpose();
    }
    if (imc>0 && imc%500==0){
      auto t2 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> deltat = t2-t1;
      std::cout<<imc<<"th mcmc steps, time elapsed: "<< deltat.count()<<std::endl;
      t1 = t2;
    }
  }

  auto t3 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> deltat_all = t3-t0;
  double RunTime = deltat_all.count();
  return List::create(
    Rcpp::Named("tau_mcmc") = tau_records,
    Rcpp::Named("sig_mcmc") = sig_records,
    Rcpp::Named("sig2_mean_mcmc") = sig2_mean_records,
    Rcpp::Named("Z_mcmc") = Z_records,
    Rcpp::Named("Accept_ratio") = n_accept / (double)(n_mcmc-n_burn),
    Rcpp::Named("L") = L,
    Rcpp::Named("D") = D,
    Rcpp::Named("RunTime") = RunTime,
    Rcpp::Named("CG_iters") = cg_iters_records,
    Rcpp::Named("CG_err") = cg_err_records
  );
}


