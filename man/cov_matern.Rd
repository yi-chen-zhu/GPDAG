% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{cov_matern}
\alias{cov_matern}
\title{Covariance matrix of Matern GP}
\usage{
cov_matern(X1, X2, nu = 3/2, tau = 1, s = 1)
}
\arguments{
\item{X1}{arma::mat, each row representing the coordiantes of a location}

\item{X2}{arma::mat, each row representing the coordiantes of a location}

\item{nu}{double, the smoothness of the Matern process}

\item{tau}{double, the time rescaling parameter of the Matern process}

\item{s}{double, the space rescaling parameter of the Matern process. Currently set as 1}
}
\value{
the computed covariance matrix.
}
\description{
This function evaluates the covariance matrix of Matern GP using the coordinates.
}
