#' Simulated effect of protein on milk yield
#'
#' Simulated milk yield in population of pedigree \code{\link[pedigreemm:pedSires]{pedSires}}.
#' Example for LMMs. Simulation based on model \code{y ~ protein + (1|herd) + (protein|sire)} with
#' \code{beta = c(6000, 200)} and \code{var_comp = 500  400  300 -200   50}
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
