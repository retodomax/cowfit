#' Simulated effect of protein on milk yield
#'
#' Simulated milk yield in population of pedigree \code{\link[pedigreemm:pedSires]{pedSires}}.
#' Example for LMMs. Simulation based on model \code{y ~ protein + (1|herd) + (protein|sire)} with
#' \code{beta = c(6000, 200)} and \code{var_comp = c(500, 400, 300, -200, 50)}.
#'
#' @docType data
#'
#' @usage data(sim_milk)
#'
#' @format An object of class \code{data.frame}
#'
#' @keywords datasets
#'
#' @examples
#' data(sim_milk)
"sim_milk"




#' Simulated Milk Fat
#'
#' Simulated milk fat in population of pedigree \code{\link[pedigreemm:pedCows]{pedCows}}.
#' Example for LMMs. Simulation based on model \code{y ~ (1|animal)} with
#' \code{beta = 850} and \code{var_comp = c(200, 60)}.
#'
#' @docType data
#'
#' @usage data(sim_fat)
#'
#' @format An object of class \code{data.frame}
#'
#' @keywords datasets
#'
#' @examples
#' data(sim_fat)
"sim_fat"




#' Simulated effect of lactation on multiple birth
#'
#' Simulated birth in population of pedigree \code{\link[pedigreemm:pedSires]{pedSires}}.
#' Example for GLMMs. Simulation based on model \code{y ~ lact + (lact|sire)} with
#' \code{beta = c(-3, 0.5)} and \code{var_comp = c(4, 3, -0.5, 1)}. \code{y = 0} stands for
#' single birth and \code{y = 1} for multiple birth.
#'
#' @docType data
#'
#' @usage data(sim_twin)
#'
#' @format An object of class \code{data.frame}
#'
#' @keywords datasets
#'
#' @examples
#' data(sim_twin)
"sim_twin"
