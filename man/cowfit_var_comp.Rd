% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cowfit_lme4.R
\name{cowfit_var_comp}
\alias{cowfit_var_comp}
\title{Variance Component Structure}
\usage{
cowfit_var_comp(formula, data, control = lmerControl())
}
\arguments{
\item{formula}{Formula which will be used in fitting function
\code{\link{cowfit_lmer}} or \code{\link{cowfit_glmer}}}

\item{data}{Data set containing all variables used in formula.
Make sure the columns have the correct data type (numeric, factor, ...)}

\item{control}{special control options of lmerControl()}
}
\value{
data frame with 4 columns
    \itemize{
      \item \code{grp} groping factor associated with variance component
      \item \code{var1} and \code{var2} names of variables associated
      with variance component (\code{var2} is \code{<NA>} unless
      it is a covariance parameter)
      \item \code{vcov} Variance or Covariance to be filled in
    }
}
\description{
The function helps to identify the necessary variance components
    for a given model. It returns a template to fill in the variance components.
    The resulting variance component vector has the correct order
    which is necessary for the argument \code{var_comp} in function
    \code{\link{cowfit_lmer}} or \code{\link{cowfit_glmer}}.
    This is especially useful for more complicated random effects models
    e.g. \code{y ~ x1 + (1 + x2 + x3:x2 | fac)}
}
\examples{
(myvar <- cowfit_var_comp(formula = y ~ (1|herd) + (protein|sire),
                          data = sim_milk))
myvar$vcov <- c(500, 400, 300, -200, 50)
cowfit_lmer(formula = y ~ (1|herd) + (protein|sire), data = sim_milk,
            pedigree = list(sire = pedSires), var_comp = myvar$vcov)
}
\seealso{
\code{\link[lme4:VarCorr]{VarCorr}}
}
\author{
Reto Zihlmann
}
