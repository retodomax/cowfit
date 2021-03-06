% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cowfit_brms.R
\name{validate_recov_matrix}
\alias{validate_recov_matrix}
\title{Internal function of \code{validate_recov_matrix}}
\usage{
validate_recov_matrix(M)
}
\description{
The orignial \code{validate_recov_matrix} expects a
    correlation matrix as input and checks it to be or symmetric.
    In this adapted version the matrix will be the cholesky factor of the
    correlation matrix and there is no check for symmetry.
}
\keyword{Internal}
