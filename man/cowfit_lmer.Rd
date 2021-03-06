% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cowfit_lme4.R
\name{cowfit_lmer}
\alias{cowfit_lmer}
\title{Fit LMMs in genetic evaluation}
\usage{
cowfit_lmer(
  formula,
  data = NULL,
  pedigree = list(),
  var_comp = NULL,
  exact_var_comp = FALSE,
  cowfit_verbose = TRUE,
  REML = TRUE,
  control = lmerControl(),
  start = NULL,
  verbose = 0L,
  subset,
  weights,
  na.action,
  offset,
  contrasts = NULL,
  devFunOnly = FALSE,
  ...
)
}
\arguments{
\item{formula}{Model formula. As in \code{\link[lme4:lmer]{lmer}}.}

\item{data}{As in \code{\link[lme4:lmer]{lmer}}.}

\item{pedigree}{a named list of
\code{\link[pedigreemm:pedigree-class]{pedigree}} objects.
The names must correspond to the names of grouping factors
for random-effects terms in the formula argument.}

\item{var_comp}{Ordered vector of variance components. Correct ordering can be
obtained with \code{\link{cowfit_var_comp}}.}

\item{exact_var_comp}{Logical. Whether numerical optimization should be used
to find theta scaling factor which will lead to exact variance components
in the final output. If \code{FALSE}, only the ratio of variance components
to redsidual variance will be conserved.}

\item{cowfit_verbose}{Logical. Verbose output is generated showing which module
of the function is currently active.}

\item{REML}{As in \code{\link[lme4:lmer]{lmer}}.}

\item{control}{As in \code{\link[lme4:lmer]{lmer}}.}

\item{start}{As in \code{\link[lme4:lmer]{lmer}}.}

\item{verbose}{As in \code{\link[lme4:lmer]{lmer}}.}

\item{subset}{As in \code{\link[lme4:lmer]{lmer}}.}

\item{weights}{As in \code{\link[lme4:lmer]{lmer}}.}

\item{na.action}{As in \code{\link[lme4:lmer]{lmer}}.}

\item{offset}{As in \code{\link[lme4:lmer]{lmer}}.}

\item{contrasts}{As in \code{\link[lme4:lmer]{lmer}}.}

\item{devFunOnly}{As in \code{\link[lme4:lmer]{lmer}}.}
}
\value{
Fitted object of class \code{\link[lmercowfit-class]{lmercowfit}}
}
\description{
Fit LMMs with random effect correlated according to pedigree.
}
\examples{
(myvar <- cowfit_var_comp(formula = y ~ (1|herd) + (protein|sire),
                          data = sim_milk))
myvar$vcov <- c(500, 400, 300, -200, 50)
fit <- cowfit_lmer(formula = y ~ (1|herd) + (protein|sire),
                   data = sim_milk,
                   pedigree = list(sire = pedSires),
                   var_comp = myvar$vcov)
as.data.frame(VarCorr(fit))
ranef(fit)
}
\seealso{
\code{\link[lme4:lmer]{lmer}}, \code{\link[pedigreemm:pedigreemm]{pedigreemm}}
}
\author{
Reto Zihlmann
}
