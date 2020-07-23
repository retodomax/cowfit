#### Import #############################################################

# We are going to import all (also internal) functions of brms.
# If we would simply use '@import brms' only the exported functions
# would be available but we need to access all of them
for(i in ls(asNamespace("brms"))){
  assign(i, utils::getFromNamespace(i, "brms"))
}



#### Exported ############################################################

# var_comp_to_stan_format() -----------------------------------------------

#' @title Transform variance component vector to list for rstan
#' @description Variance components need to be passed to rstan in the form of
#'     a list containing the standard deviations as \code{sd_i} terms and
#'     cholesky factors of correlation matrices as \code{L_i} terms.
#' @param var_comp Ordered vector of variance components. Correct ordering can be
#'     obtained with \code{\link{cowfit_var_comp}}.
#' @param formula Model formula. As in \code{\link[lme4:lmer]{lmer}}.
#' @param data As in \code{\link[lme4:lmer]{lmer}}.
#' @return List containing variance components used as input for \code{\link[rstan:sampling]{rstan::sampling}}
#' @author Reto Zihlmann
#' @seealso \code{\link[rstan:sampling]{rstan::sampling}}
#' @examples
#'     var_comp <- var_comp_to_stan_format(var_comp = c(500, 400, 300, -200, 50),
#'                                         formula = y ~ protein + (1|herd) + (protein|sire),
#'                                         data = sim_milk)
#' @export
var_comp_to_stan_format <- function(var_comp, formula, data){
  mc <- match.call()
  mc[[1]] <- as.name("lFormula")
  mc$var_comp <- NULL
  mc$formula <- update(formula, y ~.)
  data$y <- 0
  mc$data <- data
  control <-  lme4::lmerControl()
  control$checkControl$check.nobs.vs.nlev <- "ignore"
  control$checkControl$check.nobs.vs.nRE <- "ignore"
  mc[["control"]] <- control
  lmf <- eval(mc, parent.frame(1L))   # calls lFormula() with all arguments passed to cowfit_var_comp()
  cnms <- lmf$reTrms$cnms           # extract the list of all random effect Terms

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
    l_sd <- ncol(Sigma_i)
    sd_i <- list(as.array(sqrt(diag(Sigma_i))), as.matrix(t(chol(cov2cor(Sigma_i)))))
    names(sd_i) <- c(paste0("sd_", i), paste0("L_", i))
    if (l_sd < 2){
      sd_i <- sd_i[1]
    }
    out[[i]] <- sd_i
  }
  unlist(out, recursive = FALSE)
}


# cowfit_brm() ------------------------------------------------------------

#' @title Fit bayesian GLMMs in genetic evaluation
#' @description Fit Bayesian GLMMs with random effect correlated according to pedigree
#'     using Stan for full Bayesian inference.
#' @param var_comp Ordered vector of variance components. Correct ordering can be
#'     obtained with \code{\link{cowfit_var_comp}}.
#' @param formula As in \code{\link[lme4:glmer]{glmer}}.
#' @param data As in \code{\link[lme4:glmer]{glmer}}.
#' @param family As in \code{\link[lme4:glmer]{glmer}}.
#' @param pedigree a named list of
#'     \code{\link[pedigreemm:pedigree-class]{pedigree}} objects.
#'     The names must correspond to the names of grouping factors
#'     for random-effects terms in the formula argument.
#' @param var_comp Ordered vector of variance components. Correct ordering can be
#'     obtained with \code{\link{cowfit_var_comp}}.
#' @param prior As in \code{\link[brms:brm]{brm}}.
#' @param autocor As in \code{\link[brms:brm]{brm}}.
#' @param data2 As in \code{\link[brms:brm]{brm}}.
#' @param cov_ranef As in \code{\link[brms:brm]{brm}}.
#' @param sample_prior As in \code{\link[brms:brm]{brm}}.
#' @param sparse As in \code{\link[brms:brm]{brm}}.
#' @param knots As in \code{\link[brms:brm]{brm}}.
#' @param stanvars As in \code{\link[brms:brm]{brm}}.
#' @param stan_funs As in \code{\link[brms:brm]{brm}}.
#' @param fit As in \code{\link[brms:brm]{brm}}.
#' @param save_ranef As in \code{\link[brms:brm]{brm}}.
#' @param save_mevars As in \code{\link[brms:brm]{brm}}.
#' @param save_all_pars As in \code{\link[brms:brm]{brm}}.
#' @param chains As in \code{\link[brms:brm]{brm}}.
#' @param iter As in \code{\link[brms:brm]{brm}}.
#' @param warmup As in \code{\link[brms:brm]{brm}}.
#' @param thin As in \code{\link[brms:brm]{brm}}.
#' @param cores As in \code{\link[brms:brm]{brm}}.
#' @param control As in \code{\link[brms:brm]{brm}}.
#' @param algorithm As in \code{\link[brms:brm]{brm}}.
#' @param backend As in \code{\link[brms:brm]{brm}}.
#' @param future As in \code{\link[brms:brm]{brm}}.
#' @param silent As in \code{\link[brms:brm]{brm}}.
#' @param seed As in \code{\link[brms:brm]{brm}}.
#' @param save_model As in \code{\link[brms:brm]{brm}}.
#' @param stan_model_args As in \code{\link[brms:brm]{brm}}.
#' @param file As in \code{\link[brms:brm]{brm}}.
#' @param empty As in \code{\link[brms:brm]{brm}}.
#' @param rename As in \code{\link[brms:brm]{brm}}.
#' @return Fitted object of class \code{\link[brms:brmsfit-class]{brmsfit}}
#' @author Reto Zihlmann
#' @seealso \code{\link[brms:brm]{brm}}, \code{\link[cowfit_glmer]{cowfit_glmer}}
#' @examples
#' (myvar <- cowfit_var_comp(formula = y ~ (1|herd) + (protein|sire),
#'                           data = sim_milk))
#' myvar$vcov <- c(500, 400, 300, -200, 50)
#' fit <- cowfit_brm(y ~ protein + (1|herd) + (protein|sire),
#'                   data = sim_milk,
#'                   pedigree = list(sire = pedSires),
#'                   var_comp = c(500, 400, 300, -200, 50),
#'                   chains = 1)
#' @export
cowfit_brm <- function(formula, data, family = gaussian(),
                       pedigree = list(), var_comp = NULL,
                       prior = NULL, autocor = NULL,
                       data2 = NULL, cov_ranef = NULL, sample_prior = "no", sparse = NULL,
                       knots = NULL, stanvars = NULL, stan_funs = NULL, fit = NA,
                       save_ranef = TRUE, save_mevars = FALSE, save_all_pars = FALSE,
                       inits = "random", chains = 4, iter = 2000, warmup = floor(iter/2),
                       thin = 1, cores = getOption("mc.cores", 1L), control = NULL,
                       algorithm = getOption("stan_algorithm", "sampling"),
                       backend = getOption("stan_backend", "rstan"),
                       future = getOption("future", FALSE), silent = TRUE,
                       seed = NA, save_model = NULL, stan_model_args = list(), file = NULL,
                       empty = FALSE, rename = TRUE, ...) {
  mc <- mcout <- match.call()
  stopifnot(is.list(pedigree), length(names(pedigree)) == length(pedigree),
            all(sapply(pedigree, is, class2 = "pedigree")))
  if (!is.null(var_comp)){
    var_comp <- var_comp_to_stan_format(var_comp = var_comp, formula = formula, data = data)
  }
  get_relfactor <- function(name, pedigree, data) {
    mylabs <- unique(data[[name]])
    t(relfactor(ped = pedigree[[name]], labs = mylabs))
  }
  all_relfactor <- lapply(names(pedigree), get_relfactor, pedigree = pedigree, data = data)
  names(all_relfactor) <- paste0(names(pedigree), "_L")
  sub_arg <- as.list(parse(text = paste0("gr(", names(pedigree), ", cov = ", names(all_relfactor), ")")))
  names(sub_arg) <- names(pedigree)
  formula <- do.call("substitute", list(formula, sub_arg))
  data2 <- c(data2, all_relfactor)
  mc[[1]] <- quote(adapted_brm)
  mc$pedigree <- NULL
  mc$formula <- formula
  mc$data2 <- data2
  mc$var_comp <- var_comp
  fit <- eval(mc, parent.frame(1L))
  fit
}


# format_ranef() ----------------------------------------------------------

#' @title Random effects in \code{lme4} format
#' @description Returns the random effect of a fitted model in the \code{lme4}
#'     format.
#' @param object An object of a class of fitted models with random effects, including classes
#'     \code{\link[lme4:merMod-class]{merMod}}, \code{\link{cowfit}} and
#'     \code{\link[brms:brmsfit-class]{brmsfit}}.
#' @return List of data frames with the random effects of each term
#' @author Reto Zihlmann
#' @seealso \code{\link[lme4:merMod-class]{merMod}}, \code{\link{cowfit}}, \code{\link[brms:brmsfit-class]{brmsfit}}
#' @examples
#'     fit <- cowfit_lmer(y ~ (1|sire), data = sim_milk,
#'                        pedigree = list(sire = pedSires))
#'     format_ranef(fit)
#' @export
format_ranef <- function(object){
  if (class(object) == "brmsfit"){
    format_brms_ranef <- function(term_array){
      term_df <- as.data.frame(term_array[,1,,drop = FALSE])
      mynames <- dimnames(term_array)[[3]]
      colnames(term_df) <- gsub("Intercept", "(Intercept)", mynames)
      term_df
    }
    ranef_array <- ranef(object)
    lapply(ranef_array, format_brms_ranef)
  } else {
    ranef(object)
  }
}


# test_cowfit_brm() -------------------------------------------------------

#' @title Test fitting function
#' @description Test \code{\link{cowfit_brm}} for speed and accuracy of predicted
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
#' out <- test_cowfit_brm(formula = y ~ lact + (lact|sire), data = sim_twin,
#'                        family = "binomial", pedigree = list(sire = pedSires),
#'                        var_comp = c(4, 3, -0.5, 1), beta = c(-3, 0.5),
#'                        given_var_comp = TRUE)
#' @export
test_cowfit_brm <- function(formula = ~ lact + (1|herd) + (lact|sire), data = milk,
                            pedigree = list(sire = pedCowsR),
                            family = "gaussian",
                            var_comp = c(4,6,2,0.5,1), beta = c(100, 2),
                            given_var_comp = TRUE, return_all = TRUE, ...){
  mysim <- sim_glmer(formula = formula, data = data, pedigree = pedigree, family = family,
                     var_comp = var_comp, beta = beta, return_ranef = TRUE)
  if(!given_var_comp){
    var_comp <- NULL
  }
  ti <- system.time({
    fit <- cowfit_brm(update(formula, y ~ .), data = mysim$data, family = family, pedigree = pedigree,
                      var_comp = var_comp, ...)
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


#### Internal ############################################################

## brm() --> validate_data2() --> validate_recov_matrix()
##                                  xxxxxxxxxxxxxxxx
##                                  Not checking for symmetric etc. (no cov matrix but Cholesky)

## brm() --> .make_stancode() --> stan_re() --> .stan_re()
##                                              xxxxx
##                                              Assign sd_ and L_ terms to data (instead of parameters)
##                                              -- needs var_comp

## brm() --> .make_standata() --> data_gr_global()
##                                  xxxxxxxxxxx
##                                  1) Not t(chol(A)) but directly use it as Cholesky
##                                  2) Add sd's and L's in data
##                                  --- Needs var_comp


# validate_recov_matrix() -------------------------------------------------

#' @title Internal function of \code{validate_recov_matrix}
#' @description The orignial \code{validate_recov_matrix} expects a
#'     correlation matrix as input and checks it to be or symmetric.
#'     In this adapted version the matrix will be the cholesky factor of the
#'     correlation matrix and there is no check for symmetry.
#' @keywords Internal
validate_recov_matrix <- function (M) {
  M <- as.matrix(M)
  found_levels <- rownames(M)
  if (is.null(found_levels)) {
    found_levels <- colnames(M)
  }
  if (is.null(found_levels)) {
    stop2("Row or column names are required for within-group covariance matrices.")
  }
  rownames(M) <- colnames(M) <- found_levels
  M
}


# validate_data2() --------------------------------------------------------

# No changes necessary except for environmental property
environment(validate_data2) <- asNamespace("cowfit")


# data_gr_global() --------------------------------------------------------

#' @title Adapted \code{data_gr_global}
#' @description Changes the orignial \code{data_gr_global} at two places.
#'     First, the function expects a cholesky factor and not a correlation matrix
#'     in \code{data2} i.e. no cholesky decomposition needs to be done.
#'     Second, the function adds the list of variance components to
#'     the output data.
#' @keywords Internal
data_gr_global <- function (ranef, data2, var_comp = NULL) {
  out <- list()
  for (id in unique(ranef$id)) {
    tmp <- list()
    id_ranef <- subset2(ranef, id = id)
    nranef <- nrow(id_ranef)
    group <- id_ranef$group[1]
    levels <- attr(ranef, "levels")[[group]]
    tmp$N <- length(levels)
    tmp$M <- nranef
    tmp$NC <- as.integer(nranef * (nranef - 1)/2)
    if (nzchar(id_ranef$by[1])) {
      stopifnot(!nzchar(id_ranef$type[1]))
      bylevels <- id_ranef$bylevels[[1]]
      Jby <- match(attr(levels, "by"), bylevels)
      tmp$Nby <- length(bylevels)
      tmp$Jby <- as.array(Jby)
    }
    cov <- id_ranef$cov[1]
    if (nzchar(cov)) {
      ####### Changed to allow direct Cholesky (without Covariance matrix)
      Lcov <- as.matrix(data2[[cov]])
      found_levels <- rownames(Lcov)
      found <- levels %in% found_levels
      if (any(!found)) {
        stop2("Levels of the within-group covariance matrix for '",
              group, "' do not match names of the grouping levels.")
      }
      Lcov <- Lcov[levels, levels, drop = FALSE]
      tmp$Lcov <- Lcov
      #################### end of change
    }
    names(tmp) <- paste0(names(tmp), "_", id)
    c(out) <- tmp
  }
  c(out, var_comp)
}


# .make_standata() --------------------------------------------------------

#' @title Adapted \code{.make_standata()}
#' @description Changes original \code{.make_standata()} such that it takes argument
#'     \code{var_comp} and passes it to function \code{data_gr_global}.
#' @keywords Internal
.make_standata <- function (bterms, data, prior, stanvars, data2, check_response = TRUE,
                               only_response = FALSE, internal = FALSE, basis = NULL, var_comp = NULL, ...) {
  check_response <- as_one_logical(check_response)
  only_response <- as_one_logical(only_response)
  internal <- as_one_logical(internal)
  data <- order_data(data, bterms = bterms)
  out <- data_response(bterms, data, check_response = check_response,
                       internal = internal, basis = basis)
  if (!only_response) {
    ranef <- tidy_ranef(bterms, data, old_levels = basis$levels)
    meef <- tidy_meef(bterms, data, old_levels = basis$levels)
    c(out) <- data_predictor(bterms, data = data, prior = prior,
                             data2 = data2, ranef = ranef, basis = basis)
    ####################################
    c(out) <- data_gr_global(ranef, data2 = data2, var_comp = var_comp)
    ####################################
    c(out) <- data_Xme(meef, data = data)
  }
  out$prior_only <- as.integer(is_equal(get_sample_prior(prior),
                                        "only"))
  if (is.stanvars(stanvars)) {
    stanvars <- subset_stanvars(stanvars, block = "data")
    inv_names <- intersect(names(stanvars), names(out))
    if (length(inv_names)) {
      stop2("Cannot overwrite existing variables: ", collapse_comma(inv_names))
    }
    out[names(stanvars)] <- lapply(stanvars, "[[", "sdata")
  }
  if (internal) {
    attr(out, "old_order") <- attr(data, "old_order")
    attr(out, "levels") <- get_levels(tidy_meef(bterms, data),
                                      tidy_ranef(bterms, data))
  }
  structure(out, class = "standata")
}


# .stan_re() --------------------------------------------------------------

#' @title Adapted \code{.stan_re}
#' @description Changes to original \code{.stan_re}
#'     \itemize{
#'         \item Takes argument var_comp
#'         \item Puts sd_i and L_i terms to data instead of parameters if \code{!is.null(var_comp)}
#'         \item Removes prior terms for sd_i and L_i if \code{!is.null(var_comp)}
#'     }
#' @keywords Internal
.stan_re <- function (id, ranef, prior, var_comp = NULL) {
  out <- list()
  r <- subset2(ranef, id = id)
  has_cov <- nzchar(r$cov[1])
  has_by <- nzchar(r$by[[1]])
  Nby <- seq_along(r$bylevels[[1]])
  ng <- seq_along(r$gcall[[1]]$groups)
  px <- check_prefix(r)
  uresp <- usc(unique(px$resp))
  idp <- paste0(r$id, usc(combine_prefix(px)))
  str_add(out$data) <- glue("  // data for group-level effects of ID {id}\n",
                            "  int<lower=1> N_{id};  // number of grouping levels\n",
                            "  int<lower=1> M_{id};  // number of coefficients per level\n")
  if (r$gtype[1] == "mm") {
    for (res in uresp) {
      str_add(out$data) <- cglue("  int<lower=1> J_{id}{res}_{ng}[N{res}];",
                                 "  // grouping indicator per observation\n",
                                 "  real W_{id}{res}_{ng}[N{res}];", "  // multi-membership weights\n")
    }
  }
  else {
    str_add(out$data) <- cglue("  int<lower=1> J_{id}{uresp}[N{uresp}];",
                               "  // grouping indicator per observation\n")
  }
  if (has_by) {
    str_add(out$data) <- glue("  int<lower=1> Nby_{id};  // number of by-factor levels\n",
                              "  int<lower=1> Jby_{id}[N_{id}];", "  // by-factor indicator per observation\n")
  }
  if (has_cov) {
    str_add(out$data) <- glue("  matrix[N_{id}, N_{id}] Lcov_{id};",
                              "  // cholesky factor of known covariance matrix\n")
  }
  J <- seq_rows(r)
  reqZ <- !r$type %in% "sp"
  if (any(reqZ)) {
    str_add(out$data) <- "  // group-level predictor values\n"
    if (r$gtype[1] == "mm") {
      for (i in which(reqZ)) {
        str_add(out$data) <- cglue("  vector[N{usc(r$resp[i])}] Z_{idp[i]}_{r$cn[i]}_{ng};\n")
      }
    }
    else {
      str_add(out$data) <- cglue("  vector[N{usc(r$resp[reqZ])}] Z_{idp[reqZ]}_{r$cn[reqZ]};\n")
    }
  }
  if (has_by) {
    str_add_list(out) <- stan_prior(prior, class = "sd",
                                    group = r$group[1], coef = r$coef, type = glue("matrix<lower=0>[M_{id}, Nby_{id}]"),
                                    coef_type = glue("row_vector<lower=0>[Nby_{id}]"),
                                    suffix = glue("_{id}"), px = px, broadcast = "matrix",
                                    comment = "group-level standard deviations")
  }
  else {
    ############################# Here we change:

    mypri <- stan_prior(prior, class = "sd",
                        group = r$group[1], coef = r$coef, type = glue("vector<lower=0>[M_{id}]"),
                        coef_type = "real<lower=0>", suffix = glue("_{id}"),
                        px = px, comment = "group-level standard deviations")
    if(!is.null(var_comp)){
      mypri$data <- mypri$par
      mypri$par <- NULL
      mypri$prior <- NULL
    }
    str_add_list(out) <- mypri

    ###############################################
  }
  dfm <- ""
  tr <- get_dist_groups(r, "student")
  if (nrow(r) > 1L && r$cor[1]) {
    str_add(out$data) <- glue("  int<lower=1> NC_{id};  // number of group-level correlations\n")
    str_add(out$par) <- glue("  matrix[M_{id}, N_{id}] z_{id};",
                             "  // standardized group-level effects\n")
    str_add(out$prior) <- glue("  target += std_normal_lpdf(to_vector(z_{id}));\n")
    if (has_rows(tr)) {
      dfm <- glue("rep_matrix(dfm_{tr$ggn[1]}, M_{id}) .* ")
    }
    if (has_by) {
      if (has_cov) {
        stop2("Cannot combine 'by' variables with customized covariance ",
              "matrices when fitting multiple group-level effects.")
      }
      str_add_list(out) <- stan_prior(prior, class = "L",
                                      group = r$group[1], coef = Nby, type = glue("cholesky_factor_corr[M_{id}]"),
                                      coef_type = glue("cholesky_factor_corr[M_{id}]"),
                                      suffix = glue("_{id}"), dim = glue("[Nby_{id}]"),
                                      comment = "cholesky factor of correlation matrix")
      str_add(out$tpar_def) <- glue("  matrix[N_{id}, M_{id}] r_{id};  // actual group-level effects\n")
      str_add(out$tpar_comp) <- glue("  // compute actual group-level effects\n",
                                     "  r_{id} = {dfm}scale_r_cor_by(z_{id}, sd_{id}, L_{id}, Jby_{id});\n")
      str_add(out$gen_def) <- cglue("  // compute group-level correlations\n",
                                    "  corr_matrix[M_{id}] Cor_{id}_{Nby}", " = multiply_lower_tri_self_transpose(L_{id}[{Nby}]);\n",
                                    "  vector<lower=-1,upper=1>[NC_{id}] cor_{id}_{Nby};\n")
      str_add(out$gen_comp) <- stan_cor_gen_comp(glue("cor_{id}_{Nby}"),
                                                 glue("M_{id}"))
    }
    else {
      ######################## Here we change L_ terms

      mypri <- stan_prior(prior, class = "L",
                          group = r$group[1], suffix = usc(id), type = glue("cholesky_factor_corr[M_{id}]"),
                          comment = "cholesky factor of correlation matrix")
      if(!is.null(var_comp)){
        mypri$data <- mypri$par
        mypri$par <- NULL
        mypri$prior <- NULL
      }
      str_add_list(out) <- mypri

      #########################
      if (has_cov) {
        rdef <- glue("as_matrix(kronecker(Lcov_{id},",
                     " diag_pre_multiply(sd_{id}, L_{id})) *", " to_vector(z_{id}), N_{id}, M_{id})")
      }
      else {
        rdef <- glue("(diag_pre_multiply(sd_{id}, L_{id}) * z_{id})'")
      }
      str_add(out$tpar_def) <- glue("  matrix[N_{id}, M_{id}] r_{id};  // actual group-level effects\n")
      str_add(out$tpar_comp) <- glue("  // compute actual group-level effects\n",
                                     "  r_{id} = {dfm}{rdef};\n")
      str_add(out$gen_def) <- glue("  // compute group-level correlations\n",
                                   "  corr_matrix[M_{id}] Cor_{id}", " = multiply_lower_tri_self_transpose(L_{id});\n",
                                   "  vector<lower=-1,upper=1>[NC_{id}] cor_{id};\n")
      str_add(out$gen_comp) <- stan_cor_gen_comp(cor = glue("cor_{id}"),
                                                 ncol = glue("M_{id}"))
    }
    str_add(out$tpar_def) <- "  // using vectors speeds up indexing in loops\n"
    str_add(out$tpar_def) <- cglue("  vector[N_{id}] r_{idp}_{r$cn};\n")
    str_add(out$tpar_comp) <- cglue("  r_{idp}_{r$cn} = r_{id}[, {J}];\n")
  }
  else {
    str_add(out$par) <- glue("  vector[N_{id}] z_{id}[M_{id}];",
                             "  // standardized group-level effects\n")
    str_add(out$prior) <- cglue("  target += std_normal_lpdf(z_{id}[{seq_rows(r)}]);\n")
    Lcov <- str_if(has_cov, glue("Lcov_{id} * "))
    if (has_rows(tr)) {
      dfm <- glue("dfm_{tr$ggn[1]} .* ")
    }
    if (has_by) {
      str_add(out$tpar_def) <- cglue("  vector[N_{id}] r_{idp}_{r$cn};  // actual group-level effects\n")
      str_add(out$tpar_comp) <- cglue("  r_{idp}_{r$cn} = {dfm}(sd_{id}[{J}, Jby_{id}]' .* ({Lcov}z_{id}[{J}]));\n")
    }
    else {
      str_add(out$tpar_def) <- cglue("  vector[N_{id}] r_{idp}_{r$cn};  // actual group-level effects\n")
      str_add(out$tpar_comp) <- cglue("  r_{idp}_{r$cn} = {dfm}(sd_{id}[{J}] * ({Lcov}z_{id}[{J}]));\n")
    }
  }
  out
}


# stan_re() ---------------------------------------------------------------

#' @title Adapted \code{stan_re}
#' @description Changes original \code{stan_re()} such that it takes argument
#'     \code{var_comp} and passes it to function \code{.stan_re()}.
#' @keywords Internal
stan_re <- function (ranef, prior, var_comp = NULL, ...) {
  IDs <- unique(ranef$id)
  out <- list()
  tranef <- get_dist_groups(ranef, "student")
  if (has_rows(tranef)) {
    str_add(out$par) <- "  // parameters for student-t distributed group-level effects\n"
    for (i in seq_rows(tranef)) {
      g <- usc(tranef$ggn[i])
      id <- tranef$id[i]
      str_add_list(out) <- stan_prior(prior, class = "df",
                                      group = tranef$group[i], type = "real<lower=1>",
                                      suffix = g)
      str_add(out$par) <- glue("  vector<lower=0>[N_{id}] udf{g};\n")
      str_add(out$prior) <- glue("  target += inv_chi_square_lpdf(udf{g} | df{g});\n")
      str_add(out$tpar_def) <- glue("  vector[N_{id}] dfm{g};\n")
      str_add(out$tpar_comp) <- glue("  dfm{g} = sqrt(df{g} * udf{g});\n")
    }
  }
  tmp <- lapply(IDs, .stan_re, ranef = ranef, prior = prior, var_comp = var_comp,
                ...)
  out <- collapse_lists(ls = c(list(out), tmp))
  out
}


# .make_stancode() --------------------------------------------------------

#' @title Adapted \code{.make_stancode()}
#' @description Changes original \code{.make_stancode()} such that it takes argument
#'     \code{var_comp} and passes it to function \code{stan_re()}.
#' @keywords Internal
.make_stancode <- function (bterms, data, prior, stanvars, parse = getOption("parse_stancode",
                                                                                FALSE), backend = getOption("stan_backend", "rstan"), silent = TRUE,
                               save_model = NULL, var_comp = NULL, ...) {
  parse <- as_one_logical(parse)
  backend <- match.arg(backend, backend_choices())
  silent <- as_one_logical(silent)
  ranef <- tidy_ranef(bterms, data = data)
  meef <- tidy_meef(bterms, data = data)
  scode_predictor <- stan_predictor(bterms, data = data, prior = prior,
                                    ranef = ranef, meef = meef, stanvars = stanvars)
  scode_ranef <- stan_re(ranef, prior = prior, var_comp = var_comp)
  scode_llh <- stan_llh(bterms, data = data)
  scode_global_defs <- stan_global_defs(bterms, prior = prior,
                                        ranef = ranef)
  scode_Xme <- stan_Xme(meef, prior = prior)
  scode_prior <- paste0(scode_predictor$prior, scode_ranef$prior,
                        scode_Xme$prior, stan_unchecked_prior(prior))
  scode_functions <- paste0("// generated with brms ", utils::packageVersion("brms"),
                            "\n", "functions {\n", scode_global_defs$fun, collapse_stanvars(stanvars,
                                                                                            "functions"), "}\n")
  scode_data <- paste0("data {\n", scode_predictor$data, scode_ranef$data,
                       scode_Xme$data, "  int prior_only;  // should the likelihood be ignored?\n",
                       collapse_stanvars(stanvars, "data"), "}\n")
  scode_transformed_data <- paste0("transformed data {\n",
                                   scode_global_defs$tdata_def, scode_predictor$tdata_def,
                                   collapse_stanvars(stanvars, "tdata", "start"), scode_predictor$tdata_comp,
                                   collapse_stanvars(stanvars, "tdata", "end"), "}\n")
  scode_parameters <- paste0(scode_predictor$par, scode_ranef$par,
                             scode_Xme$par)
  scode_rngprior <- stan_rngprior(prior = scode_prior, par_declars = scode_parameters,
                                  gen_quantities = scode_predictor$gen_def, prior_special = attr(prior,
                                                                                                 "special"), sample_prior = get_sample_prior(prior))
  scode_parameters <- paste0("parameters {\n", scode_parameters,
                             scode_rngprior$par, collapse_stanvars(stanvars, "parameters"),
                             "}\n")
  scode_transformed_parameters <- paste0("transformed parameters {\n",
                                         scode_predictor$tpar_def, scode_ranef$tpar_def, scode_Xme$tpar_def,
                                         collapse_stanvars(stanvars, "tparameters", "start"),
                                         scode_predictor$tpar_prior, scode_ranef$tpar_prior, scode_Xme$tpar_prior,
                                         scode_predictor$tpar_comp, scode_ranef$tpar_comp, scode_Xme$tpar_comp,
                                         collapse_stanvars(stanvars, "tparameters", "end"), "}\n")
  scode_model <- paste0("model {\n", scode_predictor$model_def,
                        collapse_stanvars(stanvars, "model", "start"), scode_predictor$model_comp_basic,
                        scode_predictor$model_comp_eta_loop, scode_predictor$model_comp_dpar_link,
                        scode_predictor$model_comp_mu_link, scode_predictor$model_comp_dpar_trans,
                        scode_predictor$model_comp_mix, scode_predictor$model_comp_arma,
                        scode_predictor$model_comp_catjoin, scode_predictor$model_comp_mvjoin,
                        "  // priors including all constants\n", scode_prior,
                        "  // likelihood including all constants\n", "  if (!prior_only) {\n",
                        scode_llh, "  }\n", scode_rngprior$model, collapse_stanvars(stanvars,
                                                                                    "model", "end"), "}\n")
  scode_generated_quantities <- paste0("generated quantities {\n",
                                       scode_predictor$gen_def, scode_ranef$gen_def, scode_Xme$gen_def,
                                       scode_rngprior$gen_def, collapse_stanvars(stanvars, "genquant",
                                                                                 "start"), scode_predictor$gen_comp, scode_ranef$gen_comp,
                                       scode_rngprior$gen_comp, scode_Xme$gen_comp, collapse_stanvars(stanvars,
                                                                                                      "genquant", "end"), "}\n")
  scode <- paste0(scode_functions, scode_data, scode_transformed_data,
                  scode_parameters, scode_transformed_parameters, scode_model,
                  scode_generated_quantities)
  scode <- expand_include_statements(scode)
  if (parse) {
    scode <- parse_model(scode, backend, silent = silent)
  }
  if (is.character(save_model)) {
    cat(scode, file = save_model)
  }
  class(scode) <- c("character", "brmsmodel")
  scode
}


# adapted_brm() -----------------------------------------------------------

#' @title Adapted version of \code{brm()}
#' @description Similar to original \code{brms::brm()} function with some
#'     changes:
#'     \itemize{
#'         \item Takes argument \code{var_comp}
#'         \item Passes \code{var_comp} to function \code{\link{.make_stancode}}
#'         \item Passes \code{var_comp} to function \code{\link{.make_standata}}
#'     }
#' @export
adapted_brm <- function (formula, data, family = gaussian(), prior = NULL, autocor = NULL,
                    data2 = NULL, cov_ranef = NULL, sample_prior = "no", sparse = NULL,
                    knots = NULL, stanvars = NULL, stan_funs = NULL, fit = NA,
                    save_ranef = TRUE, save_mevars = FALSE, save_all_pars = FALSE,
                    inits = "random", chains = 4, iter = 2000, warmup = floor(iter/2),
                    thin = 1, cores = getOption("mc.cores", 1L), control = NULL,
                    algorithm = getOption("stan_algorithm", "sampling"), backend = getOption("stan_backend",
                                                                                             "rstan"), future = getOption("future", FALSE), silent = TRUE,
                    seed = NA, save_model = NULL, stan_model_args = list(), file = NULL,
                    empty = FALSE, rename = TRUE, var_comp = NULL, ...) {
  if (!is.null(file)) {
    x <- read_brmsfit(file)
    if (!is.null(x)) {
      return(x)
    }
  }
  dots <- list(...)
  algorithm <- match.arg(algorithm, algorithm_choices())
  backend <- match.arg(backend, backend_choices())
  silent <- as_one_logical(silent)
  iter <- as_one_numeric(iter)
  warmup <- as_one_numeric(warmup)
  thin <- as_one_numeric(thin)
  chains <- as_one_numeric(chains)
  cores <- as_one_numeric(cores)
  future <- as_one_logical(future) && chains > 0L
  seed <- as_one_numeric(seed, allow_na = TRUE)
  empty <- as_one_logical(empty)
  rename <- as_one_logical(rename)
  if (is.brmsfit(fit)) {
    x <- fit
    x$criteria <- list()
    sdata <- standata(x)
    model <- compiled_model(x)
  }
  else {
    formula <- validate_formula(formula, data = data, family = family,
                                autocor = autocor, sparse = sparse, cov_ranef = cov_ranef)
    family <- get_element(formula, "family")
    bterms <- brmsterms(formula)
    data_name <- substitute_name(data)
    data <- validate_data(data, bterms = bterms, knots = knots)
    attr(data, "data_name") <- data_name
    #################################
    data2 <- validate_data2(data2, bterms = bterms, get_data2_autocor(formula),
                               get_data2_cov_ranef(formula))
    #################################
    prior <- validate_prior(prior, bterms = bterms, data = data,
                            sample_prior = sample_prior)
    stanvars <- validate_stanvars(stanvars, stan_funs = stan_funs)
    x <- brmsfit(formula = formula, data = data, data2 = data2,
                 prior = prior, stanvars = stanvars, algorithm = algorithm,
                 backend = backend, family = family)
    x$ranef <- tidy_ranef(bterms, data = x$data)
    x$exclude <- exclude_pars(x, save_ranef = save_ranef,
                              save_mevars = save_mevars, save_all_pars = save_all_pars)
    #################################
    x$model <- .make_stancode(bterms, data = data, prior = prior,
                                 stanvars = stanvars, save_model = save_model, var_comp = var_comp)
    #################################
    sdata <- .make_standata(bterms, data = data, prior = prior,
                               data2 = data2, stanvars = stanvars, var_comp = var_comp)
    #################################
    if (empty) {
      return(x)
    }
    compile_args <- stan_model_args
    compile_args$model <- x$model
    compile_args$backend <- backend
    model <- do_call(compile_model, compile_args)
  }
  fit_args <- nlist(model, sdata, algorithm, backend, iter,
                    warmup, thin, chains, cores, inits, exclude = x$exclude,
                    control, future, seed, silent, ...)
  x$fit <- do_call(fit_model, fit_args)
  if (rename) {
    x <- rename_pars(x)
  }
  if (!is.null(file)) {
    write_brmsfit(x, file)
  }
  x
}
