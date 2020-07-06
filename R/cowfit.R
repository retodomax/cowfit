
# import ------------------------------------------------------------------

#' @import lme4 Matrix pedigreemm methods
NULL


# New Classes -------------------------------------------------------------

#' @title Fitted Pedigree-based mixed effects model
#' @description A mixed effects model fit by \code{\link{cowfit_lmer}} or \code{\link{cowfit_glmer}}.
#'     This class extends class \code{\link[lme4:merMod-class]{merMod}} and includes one
#'     additional slot, \code{TAt}.
#' @slot TAt sparse matrix of class \code{\link[Matrix:dtCMatrix-class]{dtCMatrix}},
#'     which was used to transform the
#'     model matrix \eqn{Z^T} and can be used to backtransform the random effects.
#' @seealso \code{\link[cowfit:glmercowfit-class]{glmercowfit}}, \code{\link[cowfit:lmercowfit-class]{lmercowfit}}
#' @examples showClass("cowfit")
#' @export
setClass("cowfit", representation = list(TAt = "dtCMatrix"),
         contains = "merMod")

#' @title Fitted Pedigree-based LMM
#' @description A mixed effects model fit by \code{\link{cowfit_lmer}}.
#'     The class inherits from class \code{\link[cowfit:cowfit-class]{cowfit}}.
#' @seealso \code{\link[cowfit:cowfit-class]{cowfit}}, \code{\link[cowfit:glmercowfit-class]{glmercowfit}}
#' @examples showClass("lmercowfit")
#' @export
setClass("lmercowfit", representation = list(resp="lmerResp"),
         contains = "cowfit")

#' @title Fitted Pedigree-based GLMM
#' @description A mixed effects model fit by \code{\link{cowfit_glmer}}.
#'     The class inherits from class \code{\link[cowfit:cowfit-class]{cowfit}}.
#' @seealso \code{\link[cowfit:cowfit-class]{cowfit}}, \code{\link[cowfit:lmercowfit-class]{lmercowfit}}
#' @examples showClass("lmercowfit")
#' @export
setClass("glmercowfit", representation = list(resp="glmResp"),
         contains = "cowfit")


# Logit function ----------------------------------------------------------

#' @title Logit function
#' @description Calculates Logarithm of the odds
#' @param p probability [0, 1]
#' @return Logarithm of the odds  [-inf, +inf]
#' @examples logit(0.75)
#' @export
logit <- function(p){
  log(p/(1-p))
}

#' @title Inverse Logit function
#' @description Calculates Probability from logarithm of the odds
#' @param eta logarithm of the odds [-inf, +inf]
#' @return probability [0, 1]
#' @examples inv_logit(12)
#' @export
inv_logit <- function(eta){
  1/(1+exp(-eta))
}


# print_cow ---------------------------------------------------------------

#' @title Print cow
#' @description Prints cow in console
#' @examples print_cow()
#' @export
print_cow <- function(){
  cat("    ^__^
    (oo)\\_______
    (__)\\       )\\/\\
        ||----w |
        ||     ||\n\n")
}

# cowfit_var_comp ---------------------------------------------------------
# Get vector of necessary variance components
# Passing variance components to function cowfit_lmer() requires to pass them in the right order
# The order can be checked with the following function.

#' @title Variance Component Structure
#' @description The function helps to identify the necessary variance components
#'     for a given model. It returns a template to fill in the variance components.
#'     The resulting variance component vector has the correct order
#'     which is necessary for the argument \code{var_comp} in function
#'     \code{\link{cowfit_lmer}} or \code{\link{cowfit_glmer}}.
#'     This is especially useful for more complicated random effects models
#'     e.g. \code{y ~ x1 + (1 + x2 + x3:x2 | fac)}
#' @param formula Formula which will be used in fitting function
#'     \code{\link{cowfit_lmer}} or \code{\link{cowfit_glmer}}
#' @param data Data set containing all variables used in formula.
#'     Make sure the columns have the correct data type (numeric, factor, ...)
#' @param control special control options of lmerControl()
#' @return data frame with 4 columns
#'     \itemize{
#'       \item \code{grp} groping factor associated with variance component
#'       \item \code{var1} and \code{var2} names of variables associated
#'       with variance component (\code{var2} is \code{<NA>} unless
#'       it is a covariance parameter)
#'       \item \code{vcov} Variance or Covariance to be filled in
#'     }
#' @author Reto Zihlmann
#' @seealso \code{\link[lme4:VarCorr]{VarCorr}}
#' @examples
#' (myvar <- cowfit_var_comp(formula = y ~ (1|herd) + (protein|sire),
#'                           data = sim_milk))
#' myvar$vcov <- c(500, 400, 300, -200, 50)
#' cowfit_lmer(formula = y ~ (1|herd) + (protein|sire), data = sim_milk,
#'             pedigree = list(sire = pedSires), var_comp = myvar$vcov)
#' @export
cowfit_var_comp <- function(formula, data, control = lmerControl()) {
  ## 1) get lfm$reTrms$cnms  ## MAKE it the easy way by also supplying the data...
  mc <- match.call()
  mc[[1]] <- as.name("lFormula")
  mc$formula <- update(formula, y ~.)
  data$y <- 0
  mc$data <- data
  control$checkControl$check.nobs.vs.nlev <- "ignore"
  control$checkControl$check.nobs.vs.nRE <- "ignore"
  mc[["control"]] <- control
  lmf <- eval(mc, parent.frame())   # calls lFormula() with all arguments passed to cowfit_var_comp()
  cnms <- lmf$reTrms$cnms           # extract the list of all random effect Terms

  ## 2) Loop trough cnms and do procedure inspired by getAnywhere(as.data.frame.VarCorr.merMod)
  cnms_names <- names(cnms)
  out <- vector("list", length = length(cnms_names))
  for(i in seq_along(cnms_names)){
    nm <- cnms[[i]]
    n <- length(nm)
    v <- diag(n)
    lt.v <- lower.tri(v, diag = FALSE)
    var1 = nm[c(seq(n), col(v)[lt.v])]
    var2 = c(rep(NA, n), nm[row(v)[lt.v]])
    out[[i]] <- data.frame(grp = cnms_names[i], var1 = var1, var2 = var2, vcov = NA)
  }

  ## 3) output rbind out
  out <- do.call(rbind, out)
  rbind(out, data.frame(grp = "Residual", var1 = NA, var2 = NA, vcov = NA))
}





# var_to_theta ------------------------------------------------------------
# Transform var_comp vector to theta
# sigma --> Sigma --> Lambda --> theta
# sigma: variance components
# Sigma: Variance covariance matrix
# Lambda: Cholesky decomposition of Sigma
# theta: non-zero elements of Lambda
# The calculation is done for each term seperately

#' @title Transform var_comp vector to theta
#' @description Theta is the vector which is used for numerical optimization in \code{lme4}.
#'     Theta contains the non-zero elements of the left Cholesky factor of Sigma.
#'     Sigma is the variance covarince matrix of the random effects.
#' @param var_comp Ordered vector of variance components. Correct ordering can be
#'     obtained with \code{\link{cowfit_var_comp}}.
#' @param cnms List with column names of raw random effects model matrix \eqn{X_i}.
#'     Usually obtained from output of \code{\link[lme4:lFormula]{lFormula}}.
#' @return \code{theta}
#' @author Reto Zihlmann
#' @seealso \code{\link[lme4:lFormula]{lFormula}}
#' @examples
#' lmod <- lFormula(y ~ (1|herd) + (protein|sire), data = sim_milk)
#' var_to_theta(var_comp = c(500, 400, 300, -200, 50),
#'              cnms = lmod$reTrms$cnms)
#' @export
var_to_theta <- function(var_comp, cnms){
  sigma_order <- function(nc){
    M <- matrix(0,nc,nc)
    diag(M) <- 1:nc
    M[lower.tri(M)] <- (nc+1):(nc*(nc+1)/2)
    M[lower.tri(M, TRUE)]
  }
  nc <- lengths(cnms)
  ncseq <- seq_along(nc)
  lt <- split(var_comp[1:(length(var_comp)-1)], rep.int(ncseq, (nc * (nc + 1))/2)) # split the variance components according to the different terms

  out <- vector("list", length = length(nc))  # Preallocation of theta_i (i: index for term)
  for(i in seq_along(nc)){
    rowIndices <- rep(1:nc[i], 1:nc[i])
    colIndices <- sequence(1:nc[i])
    template <- sparseMatrix(rowIndices, colIndices,  # make a matrix which has correct dimensions
                             x = 1)
    template@x <- as.double(lt[[i]][sigma_order(nc[i])])     # fill in sigma
    Sigma_i <- Matrix::forceSymmetric(template, uplo = "L")  # make symmetric => Sigma
    chol_sigma <- tryCatch({chol(Sigma_i)},
                           error = function(cond){stop("var_comp does not lead to positive definite matrix. Use better var_comp values.")})
    Lambdat <- t(chol_sigma/sqrt(var_comp[length(var_comp)]))            # get left Cholesky of Sigma divided by sigma_e
    out[[i]] <- Lambdat@x                                    # extract non-zero elements
  }
  unlist(out)    ## combine all theta_i to theta
}


# get_TAt() ---------------------------------------------------------------

#' @title Get transformation factor
#' @description Transformation factor is used to transform the model
#'     matrix \eqn{Z^T} in order to make the random effects independent
#'     between animals.
#' @param lmod Output of \code{\link[lme4:lFormula]{lFormula}}
#' @param pedigree a named list of
#'     \code{\link[pedigreemm:pedigree-class]{pedigree}} objects.
#'     The names must correspond to the names of grouping factors
#'     for random-effects terms in the formula argument.
#' @return Transformation factor
#' @author Reto Zihlmann
#' @seealso \code{\link[lme4:lFormula]{lFormula}}, \code{\link[pedigreemm:pedigree-class]{pedigree}}
#' @examples
#' lmod <- lFormula(y ~ (1|herd) + (protein|sire),
#'                  data = sim_milk)
#' get_TAt(lmod = lmod, pedigree = list(sire = pedSires))
#' @export
get_TAt <- function(lmod, pedigree) {
  pnms <- names(pedigree)        # Names of pedigree factors
  fl <- lmod$reTrms$flist         # Factor list
  asgn <- attr(fl, "assign")     # Which factor corresponds to which (reordered) term
  TAt_list <- vector("list", length = length(asgn))  # output list
  for(i in seq_along(asgn)){
    fac_name <- names(fl)[asgn[i]]
    p_i <- length(lmod$reTrms$cnms[[i]])                   # dim of Sigma_i
    l_i <- length(levels(lmod$reTrms$flist[[fac_name]]))   # dim of L_Ai
    on_list <- fac_name %in% pnms                         # Check if factor i is on pedigree list
    if (on_list) {
      # get fac_levels (maybe this could be simplified)
      Zt_i <- lmod$reTrms$Ztlist[[i]]          # get Zt matrix
      fac_levels <- rownames(Zt_i)[seq(1, length(rownames(Zt_i)), p_i)]  # skip some names (Replicated p_i times)
      Lt_Ai <- relfactor(pedigree[[fac_name]], fac_levels)
      TAt_list[[i]] <- kronecker(Lt_Ai, diag(p_i))
    } else {
      TAt_list[[i]] <- diag(p_i*l_i)
    }
  }
  TAt <- Matrix::bdiag(TAt_list)
  as(TAt, "dtCMatrix")
}


# Bolker exact var comp ---------------------------------------------------

#' @title Fit LMM with Exact Predefined Variance Components
#' @description Performs numerical optimization to find the
#'     scaling factor of \eqn{\theta} leading to a model with
#'     exact predefined variance components.
#' @param devfun Function which takes theta as input and returns the deviance.
#'     Usually obtained from \code{\link{lme4:modular}[mkLmerDevfun]}
#' @param var_comp Ordered vector of variance components. Correct ordering can be
#'     obtained with \code{\link{cowfit_var_comp}}.
#' @param lmod Output of \code{\link[lme4:lFormula]{lFormula}}
#' @param mcout Call which led to the output (\code{\link[base:match.call]{match.call}})
#' @return Fitted model of class \code{\link[lme4:merMod-class]{merMod}}
#' @author Reto Zihlmann
#' @seealso \code{\link{lme4:modular}[mkLmerDevfun]}, \code{\link{cowfit_var_comp}},
#'     \code{\link[lme4:lFormula]{lFormula}}, \code{\link[lme4:merMod-class]{merMod}}
#' @examples
#' lmod <- lFormula(formula = Reaction ~ Days + (Days | Subject),
#'                  data = sleepstudy)
#' devfun <- do.call(mkLmerDevfun, lmod)
#' (myvar <- cowfit_var_comp(formula = Reaction ~ Days + (Days | Subject),
#'                           data = sleepstudy))
#' myvar$vcov <- c(600, 30, 10, 600)
#' fit <- Bolker_exact_var_comp(devfun = devfun,
#'                              var_comp = myvar$vcov,
#'                              lmod = lmod)
#' as.data.frame(VarCorr(fit))
#' @export
Bolker_exact_var_comp <- function(devfun, var_comp, lmod, mcout = quote(Bolker())) {
  manual_theta <- var_to_theta(var_comp = var_comp, cnms = lmod$reTrms$cnms)  # this theta we assume based on previous variance comp
  buildMM <- function(theta){   # build model based on given theta
    ff <- devfun(theta)
    opt <- list(par=theta, fval = ff, conv = 0)
    mm <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr,
                   mc = mcout)
    mm
  }
  objfun <- function(x, target = var_comp[-length(var_comp)]){    # build model based on theta and return SS of obtained varinace components vs true.
    scaled_theta <- manual_theta*x
    mm <- buildMM(scaled_theta)
    myvcov <- as.data.frame(VarCorr(mm))$vcov
    return(sum((myvcov[-length(myvcov)] - target)^2))
  }
  opt <- optim(fn=objfun, par = 1, method = "L-BFGS-B", lower = 0)
  mm <- buildMM(manual_theta*opt$par)
  mm
}


# cowfit_lmer() -----------------------------------------------------------

#' @title Fit LMMs in genetic evaluation
#' @description Fit LMMs with random effect correlated according to pedigree.
#' @param formula Model formula. As in \code{\link[lme4:lmer]{lmer}}.
#' @param data As in \code{\link[lme4:lmer]{lmer}}.
#' @param pedigree a named list of
#'     \code{\link[pedigreemm:pedigree-class]{pedigree}} objects.
#'     The names must correspond to the names of grouping factors
#'     for random-effects terms in the formula argument.
#' @param var_comp Ordered vector of variance components. Correct ordering can be
#'     obtained with \code{\link{cowfit_var_comp}}.
#' @param exact_var_comp Logical. Whether numerical optimization should be used
#'     to find theta scaling factor which will lead to exact variance components
#'     in the final output. If \code{FALSE}, only the ratio of variance components
#'     to redsidual variance will be conserved.
#' @param cowfit_verbose Logical. Verbose output is generated showing which module
#'     of the function is currently active.
#' @param REML As in \code{\link[lme4:lmer]{lmer}}.
#' @param control As in \code{\link[lme4:lmer]{lmer}}.
#' @param start As in \code{\link[lme4:lmer]{lmer}}.
#' @param verbose As in \code{\link[lme4:lmer]{lmer}}.
#' @param subset As in \code{\link[lme4:lmer]{lmer}}.
#' @param weights As in \code{\link[lme4:lmer]{lmer}}.
#' @param na.action As in \code{\link[lme4:lmer]{lmer}}.
#' @param offset As in \code{\link[lme4:lmer]{lmer}}.
#' @param contrasts As in \code{\link[lme4:lmer]{lmer}}.
#' @param devFunOnly As in \code{\link[lme4:lmer]{lmer}}.
#' @return Fitted object of class \code{\link[lmercowfit-class]{lmercowfit}}
#' @author Reto Zihlmann
#' @seealso \code{\link[lme4:lmer]{lmer}}, \code{\link[pedigreemm:pedigreemm]{pedigreemm}}
#' @examples
#' cowfit_var_comp(formula = y ~ (1|herd) + (protein|sire),
#'                 data = sim_milk)
#' (myvar <- cowfit_var_comp(formula = y ~ (1|herd) + (protein|sire),
#'                           data = sim_milk))
#' myvar$vcov <- c(500, 400, 300, -200, 50)
#' fit <- cowfit_lmer(formula = y ~ (1|herd) + (protein|sire),
#'                    data = sim_milk,
#'                    pedigree = list(sire = pedSires),
#'                    var_comp = myvar$vcov)
#' as.data.frame(VarCorr(fit))
#' ranef(fit)
#' @export
cowfit_lmer <- function(formula, data = NULL, pedigree = list(),
                        var_comp = NULL, exact_var_comp = FALSE,
                        cowfit_verbose = TRUE, REML = TRUE, control = lmerControl(),
                        start = NULL, verbose = 0L, subset, weights, na.action,
                        offset, contrasts = NULL, devFunOnly = FALSE, ...){
  mc <- mcout <- match.call()
  missCtrl <- missing(control)
  if (!missCtrl && !inherits(control, "lmerControl")) {
    if (!is.list(control))
      stop("'control' is not a list; use lmerControl()")
    warning("passing control as list is deprecated: please use lmerControl() instead",
            immediate. = TRUE)
    control <- do.call(lmerControl, control)
  }
  mc[[1]] <- quote(lme4::lFormula)

  ## remove arguments unknown to lme4
  mc$pedigree <- NULL
  mc$var_comp <- NULL
  mc$exact_var_comp <- NULL
  mc$cowfit_verbous <- NULL
  control$checkControl$check.nobs.vs.nlev <- "ignore"  # allows for animal model
  control$checkControl$check.nobs.vs.nRE <- "ignore"
  mc$control <- control
  stopifnot(is.list(pedigree), length(names(pedigree)) == length(pedigree),
            all(sapply(pedigree, is, class2 = "pedigree")))

  ## 1. Module: lFormula
  if(cowfit_verbous) {cat(as.character(Sys.time()), "\t \t Starting module 1: Model formula\n")}
  lmod <- eval(mc, parent.frame(1L))
  mcout$formula <- lmod$formula
  lmod$formula <- NULL
  if (!is.null(var_comp)){
    if ((length(lmod$reTrms$theta)+1) != length(var_comp)) {
      stop("Object var_comp has wrong length. Use cowfit_var_comp to specify the vector.")
    }
  }


  ## Change to transformed model
  TAt <- get_TAt(lmod = lmod, pedigree = pedigree)
  lmod$reTrms$Zt <- TAt %*% lmod$reTrms$Zt


  ## 2. Module: mkLmerDevfun
  if(cowfit_verbous) {cat(as.character(Sys.time()), "\t \t Starting module 2: Get deviance function\n")}
  devfun <- do.call(mkLmerDevfun, c(lmod, list(start = start,
                                               verbose = verbose, control = control)))


  ## 3. Module: optimizeLmer + mkMerMod
  if(cowfit_verbous) {cat(as.character(Sys.time()), "\t \t Starting module 3: Optimize deviance\n")}

  if (devFunOnly)
    return(devfun)
  if (identical(control$optimizer, "none"))
    stop("deprecated use of optimizer=='none'; use NULL instead")
  mm <- if (exact_var_comp){
    if (is.null(var_comp))
      stop("Argument var_comp needs to be specified if exact_var_comp is TRUE")
    Bolker_exact_var_comp(devfun = devfun, var_comp = var_comp, lmod = lmod, mcout = mcout)
  } else {
    if (!is.null(var_comp)) {
      manual_theta <- var_to_theta(var_comp = var_comp, cnms = lmod$reTrms$cnms)
      ff <- devfun(manual_theta)
      opt <- list(par = manual_theta, fval = ff, conv = 0)
      cc <- NULL
    } else {
      opt <- if (length(control$optimizer) == 0) {
        s <- getStart(start, environment(devfun)$pp)
        list(par = s, fval = devfun(s), conv = 1000, message = "no optimization")
      } else {
        optimizeLmer(devfun, optimizer = control$optimizer, restart_edge = control$restart_edge,
                     boundary.tol = control$boundary.tol, control = control$optCtrl,
                     verbose = verbose, start = start, calc.derivs = control$calc.derivs,
                     use.last.params = control$use.last.params)
      }
      cc <- checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,
                      lbound = environment(devfun)$lower)
    }
    mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr,
             mc = mcout, lme4conv = cc)
  }

  ## 4. Module: create output
  if(cowfit_verbous) {cat(as.character(Sys.time()), "\t \t Starting module 4: Return output\n\n")}
  mm <- do.call(new, list(Class = "lmercowfit", TAt = TAt, resp = mm@resp, Gp = mm@Gp,
                          frame = mm@frame, flist = mm@flist, cnms = mm@cnms,
                          lower = mm@lower, theta = mm@theta,
                          beta = mm@beta, u = mm@u, devcomp = mm@devcomp,
                          pp = mm@pp, optinfo = mm@optinfo))
  mm@call <- evalq(mcout)
  if(cowfit_verbous) {print_cow()}
  mm
}



# ranef() -----------------------------------------------------------------

#' @export
setMethod("ranef", signature(object = "cowfit"),
          function(object, condVar = TRUE, drop = FALSE,
                   whichel = names(ans), postVar = FALSE, pedigree = TRUE,...)
          {

            if (length(L <- list(...))>0) {
              warning(paste("additional arguments to ranef.merMod ignored:",
                            paste(names(L),collapse=", ")))
            }
            if (!missing(postVar) && missing(condVar)) {
              warning(sQuote("postVar")," is deprecated: please use ",
                      sQuote("condVar")," instead")
              condVar <- postVar
            }
            ans <- object@pp$b(1) ## not always == c(matrix(unlist(getME(object,"b"))))

            ##########
            ###### Only thing I change
            if (pedigree){ans <- as.vector(t(object@TAt) %*% ans)}
            ##########

            if (!is.null(object@flist)) {
              ## evaluate the list of matrices
              levs <- lapply(fl <- object@flist, levels)
              asgn <- attr(fl, "assign")
              cnms <- object@cnms
              nc <- lengths(cnms) ## number of terms
              ## nb <- nc * lengths(levs)[asgn] ## number of cond modes per term
              nb <- diff(object@Gp)
              nbseq <- rep.int(seq_along(nb), nb)
              ml <- split(ans, nbseq)
              for (i in seq_along(ml))
                ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE,
                                  dimnames = list(NULL, cnms[[i]]))
              ## create a list of data frames corresponding to factors
              ans <- lapply(seq_along(fl),
                            function(i) {
                              m <- ml[asgn == i]
                              b2 <- vapply(m,nrow,numeric(1))
                              ub2 <- unique(b2)
                              if (length(ub2)>1)
                                stop("differing numbers of b per group")
                              ## if number of sets of modes != number of levels (e.g. Gaussian process/phyloglmm),
                              ##   generate numeric sequence for names

                              rnms <- if (ub2==length(levs[[i]])) levs[[i]] else seq(ub2)
                              data.frame(do.call(cbind, m),
                                         row.names = rnms,
                                         check.names = FALSE)
                            })
              names(ans) <- names(fl)
              # process whichel
              stopifnot(is(whichel, "character"))
              whchL <- names(ans) %in% whichel
              ans <- ans[whchL]

              ### THIS section calls internal functions (which again call more internal function)
              ### To allow for it we need better method to pass B vector to ranef()
              # if (condVar) {
              #   sigsqr <- sigma(object)^2
              #   rp <- rePos$new(object)
              #   if(any(lengths(rp$terms) > 1L)) {
              #     ## use R machinery here ...
              #     vv <- arrange.condVar(object,condVar(object, scaled=TRUE))
              #   } else {
              #     vv <- .Call(merPredDcondVar, object@pp$ptr(), as.environment(rp))
              #     vv <- lapply(vv, "*", sigsqr)
              #   }
              #   for (i in names(ans)) {
              #     attr(ans[[i]], "postVar") <- vv[[i]]
              #   }
              # }
              if (drop)
                ans <- lapply(ans, function(el)
                {
                  if (ncol(el) > 1) return(el)
                  pv <- drop(attr(el, "postVar"))
                  el <- drop(as.matrix(el))
                  if (!is.null(pv))
                    attr(el, "postVar") <- pv
                  el
                })
              class(ans) <- "ranef.mer"
            }
            ans
          }## ranef.merMod
)




# cowfit_glmer() ----------------------------------------------------------

#' @title Fit GLMMs in genetic evaluation
#' @description Fit GLMMs with random effect correlated according to pedigree.
#' @param formula Model formula. As in \code{\link[lme4:glmer]{glmer}}.
#' @param data As in \code{\link[lme4:glmer]{glmer}}.
#' @param family As in \code{\link[lme4:glmer]{glmer}}.
#' @param pedigree a named list of
#'     \code{\link[pedigreemm:pedigree-class]{pedigree}} objects.
#'     The names must correspond to the names of grouping factors
#'     for random-effects terms in the formula argument.
#' @param var_comp Ordered vector of variance components. Correct ordering can be
#'     obtained with \code{\link{cowfit_var_comp}}.
#' @param cowfit_verbose Logical. Verbose output is generated showing which module
#'     of the function is currently active.
#' @param control As in \code{\link[lme4:glmer]{glmer}}.
#' @param start As in \code{\link[lme4:glmer]{glmer}}.
#' @param verbose As in \code{\link[lme4:glmer]{glmer}}.
#' @param nAGQ As in \code{\link[lme4:glmer]{glmer}}.
#' @param subset As in \code{\link[lme4:glmer]{glmer}}.
#' @param weights As in \code{\link[lme4:glmer]{glmer}}.
#' @param na.action As in \code{\link[lme4:glmer]{glmer}}.
#' @param offset As in \code{\link[lme4:glmer]{glmer}}.
#' @param contrasts As in \code{\link[lme4:glmer]{glmer}}.
#' @param mustart As in \code{\link[lme4:glmer]{glmer}}.
#' @param etastart As in \code{\link[lme4:glmer]{glmer}}.
#' @param devFunOnly As in \code{\link[lme4:glmer]{glmer}}.
#' @return Fitted object of class \code{\link[glmercowfit-class]{glmercowfit}}
#' @author Reto Zihlmann
#' @seealso \code{\link[lme4:glmer]{glmer}}, \code{\link[pedigreemm:pedigreemm]{pedigreemm}}
#' @examples
#' cowfit_var_comp(formula = y ~ lact + (lact|sire),
#'                 data = sim_twin)
#' (myvar <- cowfit_var_comp(formula = y ~ lact + (lact|sire),
#'                           data = sim_twin))
#' myvar$vcov <- c(4, 3, -0.5, 1)
#' fit <- cowfit_glmer(formula = y ~ lact + (lact|sire),
#'                     data = sim_twin, family = "binomial",
#'                     pedigree = list(sire = pedSires),
#'                     var_comp = myvar$vcov)
#' as.data.frame(VarCorr(fit))
#' ranef(fit)
#' @export
cowfit_glmer <- function(formula, data = NULL, family = gaussian,
                         pedigree = list(),
                         var_comp = NULL,
                         cowfit_verbous = TRUE,
                         control = glmerControl(),
                         start = NULL, verbose = 0L, nAGQ = 1L, subset, weights, na.action,
                         offset, contrasts = NULL, mustart, etastart, devFunOnly = FALSE) {
  if (!inherits(control, "glmerControl")) {
    if (!is.list(control))
      stop("'control' is not a list; use glmerControl()")
    if (class(control)[1] == "lmerControl") {
      warning("please use glmerControl() instead of lmerControl()",
              immediate. = TRUE)
      control <- c(control[!names(control) %in% c("checkConv",
                                                  "checkControl")], control$checkControl, control$checkConv)
      control["restart_edge"] <- NULL
    }
    else {
      msg <- "Use control=glmerControl(..) instead of passing a list"
      if (length(cl <- class(control))) {
        msg <- paste(msg, "of class", dQuote(cl[1]))
      }
      warning(msg, immediate. = TRUE)
    }
    control <- do.call(glmerControl, control)
  }
  mc <- mcout <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame(2))
  if (is.function(family))
    family <- family()
  if (isTRUE(all.equal(family, gaussian()))) {
    warning("calling cowfit_glmer() with family=gaussian (identity link) as a shortcut to cowfit_lmer() is deprecated;",
            " please call cowfit_lmer() directly")
    mc[[1]] <- quote(cowfit_lmer)
    mc["family"] <- NULL
    return(eval(mc, parent.frame()))
  }
  mc[[1]] <- quote(lme4::glFormula)

  ## remove arguments unknown to lme4
  mc$pedigree <- NULL
  mc$var_comp <- NULL
  mc$cowfit_verbous <- NULL
  control$checkControl$check.nobs.vs.nlev <- "ignore"  # allows for animal model
  control$checkControl$check.nobs.vs.nRE <- "ignore"
  stopifnot(is.list(pedigree), length(names(pedigree)) == length(pedigree),
            all(sapply(pedigree, is, class2 = "pedigree")))

  ## 1. Module: lFormula
  if(cowfit_verbous) {cat(as.character(Sys.time()), "\t \t Starting module 1: Model formula\n")}
  glmod <- eval(mc, parent.frame(1L))
  if (!is.null(var_comp)){
    if ((length(glmod$reTrms$theta)+1) != length(var_comp)) {
      stop("Object var_comp has wrong length. Use cowfit_var_comp to specify the vector.")
    }
  }
  mcout$formula <- glmod$formula
  glmod$formula <- NULL
  nAGQinit <- if (control$nAGQ0initStep)
    0L
  else 1L

  ## Change to transformed model
  TAt <- get_TAt(lmod = glmod, pedigree = pedigree)
  glmod$reTrms$Zt <- TAt %*% glmod$reTrms$Zt

  ## 2. Module: mkLmerDevfun
  if(cowfit_verbous) {cat(as.character(Sys.time()), "\t \t Starting module 2: Get deviance function\n")}
  devfun <- do.call(mkGlmerDevfun, c(glmod, list(verbose = verbose,
                                                 control = control, nAGQ = nAGQinit)))

  ## 3. Module: optimizeGlmer + mkMerMod
  if(cowfit_verbous) {cat(as.character(Sys.time()), "\t \t Starting module 3: Optimize deviance\n")}

  if (nAGQ == 0 && devFunOnly)
    return(devfun)
  if (is.list(start)) {
    start.bad <- setdiff(names(start), c("theta", "fixef"))
    if (length(start.bad) > 0) {
      stop(sprintf("bad name(s) for start vector (%s); should be %s and/or %s",
                   paste(start.bad, collapse = ", "), shQuote("theta"),
                   shQuote("fixef")), call. = FALSE)
    }
    if (!is.null(start$fixef) && nAGQ == 0)
      stop("should not specify both start$fixef and nAGQ==0")
    if (!is.null(start$theta) && !is.null(var_comp))
      stop("should not specify both start$theta and var_comp")
  }

  if (!is.null(var_comp)) {
    # initial Step
    if (control$nAGQ0initStep) {
      manual_theta <- var_to_theta(var_comp = var_comp, cnms = glmod$reTrms$cnms)
      ff_theta <- devfun(manual_theta)
      opt <- list(par = manual_theta, fval = ff_theta, conv = 0) # approx opt if nAGQ = 0
    }
    if (nAGQ > 0L) {
      # update GlmerDevfun()
      devfun <- updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)
      if (devFunOnly)
        return(devfun)

      # Make wrapper
      myrho <- environment(devfun)
      para_length <- length(myrho$pp$theta) + length(myrho$pp$beta0)
      theta_beta <- rep(0, para_length)
      theta_beta[myrho$dpars] <- manual_theta
      devfun_beta <- function(mybeta){
        theta_beta[-myrho$dpars] <- mybeta
        devfun(theta_beta)
      }

      start <- if(!is.null(start)){
        start
      } else {
        myrho$pp$delb
      }

      devfun_beta(start) ## this somehow seems necessary to update environment(devfun)
      ## otherwise results are much worse!
      ## Why is this necessary?

      opt <- optim(par = start, fn = devfun_beta, method = "Nelder-Mead")

      theta_beta[-myrho$dpars] <- opt$par
      ff <- devfun(theta_beta)
      opt <- list(par = theta_beta, fval = ff, conv = 0)
    } else {
      myrho <- environment(devfun)
      myrho$nAGQ <- nAGQ    ## seems strange but is like this passed to rho in mkMerMod
    }
  } else {
    ### Usual procedure
    if (control$nAGQ0initStep) {
      opt <- optimizeGlmer(devfun, optimizer = control$optimizer[[1]],
                           restart_edge = if (nAGQ == 0)
                             control$restart_edge
                           else FALSE, boundary.tol = if (nAGQ == 0)
                             control$boundary.tol
                           else 0, control = control$optCtrl, start = start,
                           nAGQ = 0, verbose = verbose, calc.derivs = FALSE)
    }
    if (nAGQ > 0L) {
      devfun <- updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)
      if (control$nAGQ0initStep) {
        start <- lme4:::updateStart(start, theta = opt$par)
      }
      if (devFunOnly)
        return(devfun)
      opt <- optimizeGlmer(devfun, optimizer = control$optimizer[[2]],
                           restart_edge = control$restart_edge, boundary.tol = control$boundary.tol,
                           control = control$optCtrl, start = start, nAGQ = nAGQ,
                           verbose = verbose, stage = 2, calc.derivs = control$calc.derivs,
                           use.last.params = control$use.last.params)
    }
  }

  cc <- if (!control$calc.derivs || !is.null(var_comp))
    NULL
  else {
    if (verbose > 10)
      cat("checking convergence\n")
    checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,
              lbound = environment(devfun)$lower)
  }
  mm <- mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr,
                 mc = mcout, lme4conv = cc)


  ## 4. Module: create output...
  if(cowfit_verbous) {cat(as.character(Sys.time()), "\t \t Starting module 4: Return output\n\n")}
  mm <- do.call(new, list(Class = "glmercowfit", TAt = TAt, resp = mm@resp, Gp = mm@Gp,
                          frame = mm@frame, flist = mm@flist, cnms = mm@cnms,
                          lower = mm@lower, theta = mm@theta,
                          beta = mm@beta, u = mm@u, devcomp = mm@devcomp,
                          pp = mm@pp, optinfo = mm@optinfo))
  mm@call <- evalq(mcout)
  if(cowfit_verbous) {print_cow()}
  mm
}



# ranef_sim ---------------------------------------------------------------

#' @title ranef_sim
#' @description Helper function for the output module of sim_lmer()
#' @export
ranef_sim <- function(b, flist, cnms, Gp){
  ans <- b

  if (!is.null(flist)) {
    ## evaluate the list of matrices
    levs <- lapply(fl <- flist, levels)
    asgn <- attr(fl, "assign")
    nc <- lengths(cnms) ## number of terms
    ## nb <- nc * lengths(levs)[asgn] ## number of cond modes per term
    nb <- diff(Gp)
    nbseq <- rep.int(seq_along(nb), nb)
    ml <- split(ans, nbseq)
    for (i in seq_along(ml))
      ml[[i]] <- matrix(ml[[i]], ncol = nc[i], byrow = TRUE,
                        dimnames = list(NULL, cnms[[i]]))
    ## create a list of data frames corresponding to factors
    ans <- lapply(seq_along(fl),
                  function(i) {
                    m <- ml[asgn == i]
                    b2 <- vapply(m,nrow,numeric(1))
                    ub2 <- unique(b2)
                    if (length(ub2)>1)
                      stop("differing numbers of b per group")
                    ## if number of sets of modes != number of levels (e.g. Gaussian process/phyloglmm),
                    ##   generate numeric sequence for names

                    rnms <- if (ub2==length(levs[[i]])) levs[[i]] else seq(ub2)
                    data.frame(do.call(cbind, m),
                               row.names = rnms,
                               check.names = FALSE)
                  })
    names(ans) <- names(fl)
    class(ans) <- "ranef.mer"
  }
  ans
}


# sim_lmer() --------------------------------------------------------------

#' @title Simulate LMMs in genetic evaluation
#' @description Simulate LMM with genetically correlated random effects
#' @param formula Model formula
#' @param data Data frame containing all predictors in \code{formula}
#' @param pedigree a named list of
#'     \code{\link[pedigreemm:pedigree-class]{pedigree}} objects.
#'     The names must correspond to the names of grouping factors
#'     for random-effects terms in the formula argument.
#' @param var_comp Ordered vector of variance components. Correct ordering can be
#'     obtained with \code{\link{cowfit_var_comp}}.
#' @param beta Ordered vector of fixed effects
#' @param retrun_ranef Logical. Should true random effects returned together with data
#' @return List containing simulated data (and optionally random effects)
#' @author Reto Zihlmann
#' @seealso \code{\link{cowfit_lmer}}
#' @examples
#' set.seed(1)
#' sim_lmer(formula = y ~ protein + (1|herd) + (protein|sire),
#'          data = sim_milk,
#'          pedigree = list(sire = pedSires),
#'          var_comp = c(500, 400, 300, -200, 50),
#'          beta = c(6000, 200), return_ranef = FALSE)
#' @export
sim_lmer <- function(formula = y ~ (1|herd) + (lact|id), data = milk, pedigree = list(),
                     var_comp, beta, return_ranef = TRUE){
  mc <- match.call()
  mc$formula <- update(formula, y ~.) # y should always be response
  data$y <- 0                         # pseudo data
  mc$data <- data                     # add new data to call
  mc[[1]] <- quote(lme4::lFormula)

  ## remove arguments unknown to lme4
  mc$pedigree <- NULL
  mc$var_comp <- NULL
  mc$beta <- NULL
  mc$return_ranef <- NULL
  control <- lmerControl()
  control$checkControl$check.nobs.vs.nlev <- "ignore"  # allows for animal model
  control$checkControl$check.nobs.vs.nRE <- "ignore"
  mc$control <- control
  stopifnot(is.list(pedigree), length(names(pedigree)) == length(pedigree),
            all(sapply(pedigree, is, class2 = "pedigree")))

  ## 1. Module: lFormula
  mc$formula <- update(formula, y ~ .)  # always y as response
  data$y <- 0                           # add some response vector
  mc$data <- data                       # update data
  lmod <- eval(mc, parent.frame(1L))
  if ((length(lmod$reTrms$theta)+1) != length(var_comp)) {
    stop("Object var_comp has wrong length. Use cowfit_var_comp to specify the vector.")
  }

  ## 2. Module: Transform Zt
  TAt <- get_TAt(lmod = lmod, pedigree = pedigree)
  ZtStar <- TAt %*% lmod$reTrms$Zt

  ## 3. Module: Sample and calculate response
  theta <- var_to_theta(var_comp = var_comp, cnms = lmod$reTrms$cnms)
  sig_e <- sqrt(var_comp[length(var_comp)])
  Lambdat <- lmod$reTrms$Lambdat
  Lambdat@x[] <- theta[lmod$reTrms$Lind]  ## x[] makes sure that the dimensions are conserved
  u <- rnorm(nrow(Lambdat), 0, sig_e)
  epsilon <- rnorm(nrow(data), 0, sig_e)
  bStar <- as.vector(t(Lambdat) %*% u)
  b <- as.vector(t(TAt) %*% bStar)
  X <- lmod$X
  data$y <- as.vector(X %*% beta + t(ZtStar) %*% bStar + epsilon)

  ## 4. Module: create Output
  if (return_ranef) {
    myranef <- ranef_sim(b = b, flist = lmod$reTrms$flist, cnms = lmod$reTrms$cnms, Gp = lmod$reTrms$Gp)
    list(data = data, ranef = myranef)
  } else {
    data
  }

}



# test_cowfit_lmer() ------------------------------------------------------

#' @title Test fitting function
#' @description Test \code{\link{cowfit_lmer}} for speed and accuracy of predicted
#'     random effects.
#' @param formula Model formula
#' @param data Data frame containing all predictors in \code{formula}
#' @param pedigree a named list of
#'     \code{\link[pedigreemm:pedigree-class]{pedigree}} objects.
#'     The names must correspond to the names of grouping factors
#'     for random-effects terms in the formula argument.
#' @param var_comp Ordered vector of variance components. Correct ordering can be
#'     obtained with \code{\link{cowfit_var_comp}}.
#' @param beta Ordered vector of fixed effects
#' @param given_var_comp Logical. Should true variance components be passed to
#'     fitting function \code{\link{cowfit_lmer}}?
#' @param exact_var_comp Logical. Should numerical optimization be used to find exact theta
#' @param return_all Logical. If \code{TRUE} not only fitted model but also
#'     \code{\link{system.time}} and correlation between true and fitted random effects.
#' @return list containing fitted model (and optionally \code{\link{system.time}} and
#'     correlation between true and fitted random effects.)
#' @author Reto Zihlmann
#' @seealso \code{\link{cowfit_lmer}}, \code{\link{sim_lmer}}
#' @examples
#' cowfit_var_comp(formula = ~ calvingYear + (1|herd) + (mastitis|sire),
#'                 data = mastitis)
#' out <- test_cowfit_lmer(formula = ~ calvingYear + (1|herd) + (mastitis|sire),
#'                         data = mastitis, pedigree = list(sire = pedSires),
#'                         var_comp = c(4, 5, 2, 0.5, 0.5),
#'                         beta = c(100, 2, 1, 3, 4, 3),
#'                         given_var_comp = TRUE, exact_var_comp = TRUE)
#' as.data.frame(VarCorr(out$fit))
#' @export
test_cowfit_lmer <- function(formula = ~lact + (1|herd) + (lact|sire), data = milk,
                             pedigree = list(sire = pedCowsR),
                             var_comp = c(4,6,2,0.5,1), beta = c(100, 2),
                             given_var_comp = TRUE, exact_var_comp = FALSE,
                             return_all = TRUE, ...){
  mysim <- sim_lmer(formula = formula, data = data, pedigree = pedigree,
                    var_comp = var_comp, beta = beta, return_ranef = TRUE)
  if(!given_var_comp){
    var_comp <- NULL
  }
  ti <- system.time({
    fit <- cowfit_lmer(update(formula, y ~ .), data = mysim$data, pedigree = pedigree,
                       var_comp = var_comp, exact_var_comp = exact_var_comp,
                       cowfit_verbous = TRUE, ...)
  })
  if(!return_all){
    return(fit)
  }
  true_ranef <- unlist(mysim$ranef, recursive = FALSE)
  est_ranef <- unlist(ranef(fit), recursive = FALSE)
  sp_cor <- vector("numeric", length = length(true_ranef))
  for(i in seq_along(true_ranef)){
    sp_cor[i] <- cor(true_ranef[[i]], est_ranef[[i]], method = "spearman")
  }
  names(sp_cor) <- names(true_ranef)
  list(ti = ti[3], sp_cor = sp_cor, fit = fit)
}



# sim_glmer() -------------------------------------------------------------
## important: residual variance component should alsways be 1

#' @title Simulate GLMMs in genetic evaluation
#' @description Simulate GLMM with genetically correlated random effects
#' @param formula Model formula
#' @param data Data frame containing all predictors in \code{formula}
#' @param family A GLM family, see \code{\link[stats:glm]{glm}}
#'     and \code{\link[stats:family]{family}}.
#' @param pedigree a named list of
#'     \code{\link[pedigreemm:pedigree-class]{pedigree}} objects.
#'     The names must correspond to the names of grouping factors
#'     for random-effects terms in the formula argument.
#' @param var_comp Ordered vector of variance components. Correct ordering can be
#'     obtained with \code{\link{cowfit_var_comp}}.
#' @param beta Ordered vector of fixed effects
#' @param retrun_ranef Logical. Should true random effects returned together with data
#' @return List containing simulated data (and optionally random effects)
#' @author Reto Zihlmann
#' @seealso \code{\link{cowfit_glmer}}
#' @examples
#' set.seed(1)
#' sim_glmer(formula = y ~ lact + (lact|sire), data = sim_twin,
#'           family = "binomial", pedigree = list(sire = pedSires),
#'           var_comp = c(4, 3, -0.5, 1), beta = c(-3, 0.5),
#'           return_ranef = FALSE)
#' @export
sim_glmer <- function(formula = y ~ (1|herd) + (lact|id), data = milk,
                      family = "binomial", pedigree = list(),
                      var_comp, beta, return_ranef = TRUE){
  mc <- match.call()
  mc$formula <- update(formula, y ~.) # y should always be response
  data$y <- 0                         # pseudo data
  mc$data <- data                     # add new data to call
  mc[[1]] <- quote(lme4::glFormula)

  ## remove arguments unknown to lme4
  mc$pedigree <- NULL
  mc$var_comp <- NULL
  mc$beta <- NULL
  mc$return_ranef <- NULL
  control <- lmerControl()
  control$checkControl$check.nobs.vs.nlev <- "ignore"  # allows for animal model
  control$checkControl$check.nobs.vs.nRE <- "ignore"
  mc$control <- control
  stopifnot(is.list(pedigree), length(names(pedigree)) == length(pedigree),
            all(sapply(pedigree, is, class2 = "pedigree")))

  ## 1. Module: lFormula
  mc$formula <- update(formula, y ~ .)  # always y as response
  data$y <- 0                           # add some response vector
  mc$data <- data                       # update data
  glmod <- eval(mc, parent.frame(1L))
  if ((length(glmod$reTrms$theta)+1) != length(var_comp)) {
    stop("Object var_comp has wrong length. Use cowfit_var_comp to specify the vector.")
  }

  ## 2. Module: Transform Zt
  TAt <- get_TAt(lmod = glmod, pedigree = pedigree)
  ZtStar <- TAt %*% glmod$reTrms$Zt

  ## 3. Module: Sample and calculate response
  theta <- var_to_theta(var_comp = var_comp, cnms = glmod$reTrms$cnms)
  sig_e <- sqrt(var_comp[length(var_comp)])
  Lambdat <- glmod$reTrms$Lambdat
  Lambdat@x[] <- theta[glmod$reTrms$Lind]  ## x[] makes sure that the dimensions are conserved
  u <- rnorm(nrow(Lambdat), 0, sig_e)
  epsilon <- rnorm(nrow(data), 0, sig_e)
  bStar <- as.vector(t(Lambdat) %*% u)
  b <- as.vector(t(TAt) %*% bStar)
  X <- glmod$X
  data$lin_pred <- as.vector(X %*% beta + t(ZtStar) %*% bStar + epsilon)

  ## Now we need to create y from the lin_pred
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  rsample <- switch (family$family,
                     binomial = function(n, mu){rbinom(n = n, size = 1, prob = mu)},
                     poisson = function(n, mu){rpois(n = n, lambda = mu)},
                     Gamma = function(n, mu){rgamma(n = n, shape = mu, rate = 1)})
  if (is.null(rsample))
  {stop("sim_glmer() currently only supports 'binomial', 'poisson' and 'Gamma' for family argument")}
  data$mu <- family$linkinv(data$lin_pred)
  data$y <- rsample(n = nrow(data), mu = data$mu)


  ## 4. Module: create Output
  if (return_ranef) {
    myranef <- ranef_sim(b = b, flist = glmod$reTrms$flist, cnms = glmod$reTrms$cnms, Gp = glmod$reTrms$Gp)
    list(data = data, ranef = myranef)
  } else {
    data
  }
}


# test_cowfit_glmer() -----------------------------------------------------

#' @title Test fitting function
#' @description Test \code{\link{cowfit_glmer}} for speed and accuracy of predicted
#'     random effects.
#' @param formula Model formula
#' @param data Data frame containing all predictors in \code{formula}
#' @param pedigree a named list of
#'     \code{\link[pedigreemm:pedigree-class]{pedigree}} objects.
#'     The names must correspond to the names of grouping factors
#'     for random-effects terms in the formula argument.
#' @param family A GLM family, see \code{\link[stats:glm]{glm}}
#'     and \code{\link[stats:family]{family}}.
#' @param var_comp Ordered vector of variance components. Correct ordering can be
#'     obtained with \code{\link{cowfit_var_comp}}.
#' @param beta Ordered vector of fixed effects
#' @param given_var_comp Logical. Should true variance components be passed to
#'     fitting function \code{\link{cowfit_glmer}}?
#' @param return_all Logical. If \code{TRUE} not only fitted model but also
#'     \code{\link{system.time}} and correlation between true and fitted random effects.
#' @return list containing fitted model (and optionally \code{\link{system.time}} and
#'     correlation between true and fitted random effects.)
#' @author Reto Zihlmann
#' @seealso \code{\link{cowfit_glmer}}, \code{\link{sim_glmer}}
#' @examples
#' cowfit_var_comp(formula = y ~ lact + (lact|sire), data = sim_twin)
#' out <- test_cowfit_glmer(formula = y ~ lact + (lact|sire), data = sim_twin,
#'                          family = "binomial", pedigree = list(sire = pedSires),
#'                          var_comp = c(4, 3, -0.5, 1), beta = c(-3, 0.5),
#'                          given_var_comp = TRUE)
#' as.data.frame(VarCorr(out$fit))
#' @export
test_cowfit_glmer <- function(formula = ~lact + (1|herd) + (lact|sire), data = milk,
                              pedigree = list(sire = pedCowsR),
                              family = "binomial",
                              var_comp = c(4,6,2,0.5,1), beta = c(100, 2),
                              given_var_comp = TRUE, return_all = TRUE, ...){
  mysim <- sim_glmer(formula = formula, data = data, pedigree = pedigree, family = family,
                     var_comp = var_comp, beta = beta, return_ranef = TRUE)
  if(!given_var_comp){
    var_comp <- NULL
  }
  ti <- system.time({
    fit <- cowfit_glmer(update(formula, y ~ .), data = mysim$data, pedigree = pedigree, family = family,
                        var_comp = var_comp, cowfit_verbous = TRUE, ...)
  })
  if(!return_all){
    return(fit)
  }
  true_ranef <- unlist(mysim$ranef, recursive = FALSE)
  est_ranef <- unlist(ranef(fit), recursive = FALSE)
  sp_cor <- vector("numeric", length = length(true_ranef))
  for(i in seq_along(true_ranef)){
    sp_cor[i] <- cor(true_ranef[[i]], est_ranef[[i]], method = "spearman")
  }
  names(sp_cor) <- names(true_ranef)
  list(ti = ti[3], sp_cor = sp_cor, fit = fit)
}
