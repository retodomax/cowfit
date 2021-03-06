% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{sim_fat}
\alias{sim_fat}
\title{Simulated Milk Fat}
\format{
An object of class \code{data.frame}
}
\usage{
data(sim_fat)
}
\description{
Simulated milk fat in population of pedigree \code{\link[pedigreemm:pedCows]{pedCows}}.
Example for LMMs. Simulation based on model \code{y ~ (1|animal)} with
\code{beta = 850} and \code{var_comp = c(200, 60)}.
}
\examples{
data(sim_fat)
}
\keyword{datasets}
