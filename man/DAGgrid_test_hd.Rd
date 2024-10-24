% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DAGgrid.R
\name{DAGgrid_test_hd}
\alias{DAGgrid_test_hd}
\title{Norming DAG for test data, 2-d grid algorithm}
\usage{
DAGgrid_test_hd(Xcoords, Xtest, alpha, sort_in_ord)
}
\arguments{
\item{Xcoords}{List of numeric vectors, each element of the list contains the coordinate values of training data}

\item{Xtest}{A numeric matrix, each row corresponds to the coordinates of a test location}

\item{alpha}{A numeric value, the smoothness of the Matern process, determining the cardinality of parent sets.}

\item{sort_in_ord}{A numeric vector representing the ordering of elements in Xcoords. One can simply input the "sort_in_ord" variable from the output of the function DAGgrid_per_hd}
}
\value{
A list of numeric vectors representing the test DAG. Each element of the list contains the indices of parents in DAG ordering.
}
\description{
This function constructs norming DAG on a test set using 2-d grid algorithm. Only
the training set is required to be on a grid, while test set can be arbitrary locations.
}
