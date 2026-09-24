inv <- function(x){
  1/x
}

#' Function to compare evzinb or evinb models with zinb and nb models
#'
#' @param object A fitted evzinb or evinb model object
#' @param winsorize Should winsorizing be done in the comparisons?
#' @param razorize Should razorizing (trimming) be done in the comparisons?
#' @param cutoff_value Integer: Which observation should be used as a basis for winsorizing/razorising. E.g. 10 means that everything larger than the 10th observation will be winsorized/razorised
#' @param init_theta Optional initial value for theta in the NB specification
#' @param nb_comparison Should comparison be made with a negative binomial model?
#' @param zinb_comparison Should comparisons be made with the zinb model? Not
#'   available for \code{evinb} objects (there is no zero-inflation component); it
#'   defaults to \code{FALSE} for those, with a message, and errors if set to
#'   \code{TRUE} explicitly.
#' @param poisson_comparison Should comparison be made with a Poisson model?
#'   Defaults to \code{TRUE} when \code{object} was itself fitted with
#'   \code{family = "poisson"}, \code{FALSE} otherwise (round9 E.3);
#'   user-overridable either way.
#' @param zip_comparison Should comparisons be made with a zero-inflated
#'   Poisson (ZIP) model? Same default as \code{poisson_comparison}. Not
#'   available for \code{evinb} objects, with the same
#'   defaults-to-\code{FALSE}-with-a-message / errors-if-\code{TRUE} behaviour
#'   as \code{zinb_comparison}.
#' @inheritParams evzinb
#'
#' @inheritSection evzinb Parallel processing
#' @inheritSection evzinb Reproducibility
#'
#' @return An object of class \code{evzinbcomp}: a list whose first element
#'   \code{model} is the original evzinb/evinb model (also available as
#'   \code{evzinb} for backwards compatibility), followed by the compared
#'   \code{nb} / \code{zinb} models (and their winsorized/razorized variants when
#'   requested).
#' @export
#'
#' @examples
#' \donttest{
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 5)
#' compare_models(model)
#' }
#'
#' \donttest{
#' data(hks)
#' hks_mod <- evzinb(osvAll ~ troopLag + lntpop + brv_AllLag_log + osvAllLagDum,
#'                   formula_pareto = ~ log1p(troopLag),
#'                   data = hks, n_bootstraps = 2, multicore = FALSE)
#' cmp <- compare_models(hks_mod)
#' compare_fit(cmp)
#' }
compare_models <- function(object, nb_comparison = TRUE, zinb_comparison = TRUE,
                           poisson_comparison = identical((object$family %||% evinf_family())$count, "poisson"),
                           zip_comparison = identical((object$family %||% evinf_family())$count, "poisson"),
                           winsorize = FALSE, razorize = FALSE, cutoff_value=10, init_theta=NULL, multicore = NULL, ncores=NULL){

  if(!inherits(object, c('evzinb','evinb'))){
    stop('compare_models() requires a fitted evzinb or evinb object.')
  }

  if(inherits(object, 'evinb')){
    if(missing(zinb_comparison)){
      message('compare_models(): evinb models have no zero-inflation component; setting zinb_comparison = FALSE.')
      zinb_comparison <- FALSE
    }else if(isTRUE(zinb_comparison)){
      stop('compare_models(): zinb_comparison is not available for evinb objects. There is no zero-inflation component to compare.')
    }
    if(missing(zip_comparison)){
      message('compare_models(): evinb models have no zero-inflation component; setting zip_comparison = FALSE.')
      zip_comparison <- FALSE
    }else if(isTRUE(zip_comparison)){
      stop('compare_models(): zip_comparison is not available for evinb objects. There is no zero-inflation component to compare.')
    }
  }

  # round10 0.2 (review §2): every competitor fit was unweighted, so a
  # weighted evzinb()/evinb() fit's AIC/BIC comparisons were between a
  # weighted and an unweighted likelihood. object$weights_col (NULL when
  # unweighted) already names a column of object$data$data -- and of every
  # data frame derived from it below, since winsorising/razorising only
  # mutate or subset rows -- so cmp_fit() below references it by name via
  # do.call() rather than passing the weight *vector* directly: glm()/
  # glm.nb()/zeroinfl() re-resolve `weights` against environment(formula)
  # when it is not found in `data`, not against the caller's frame, once
  # model.frame()'s NSE machinery runs inside the try() these are wrapped
  # in, so a vector held in a local variable fails to resolve there.
  wcol <- object$weights_col
  cmp_fit <- function(fn, args) {
    if (!is.null(wcol)) {
      args$weights <- as.name(wcol)
    }
    try(do.call(fn, args))
  }

  deparse_rhs <- function(f){
    if(is.null(f)) return('1')
    rhs <- if(length(f) == 3L) f[[3]] else f[[2]]
    paste(deparse(rhs), collapse = ' ')
  }

  dv_f <- all.vars(object$formulas$formula_nb)[1]
  rhs_nb <- deparse_rhs(object$formulas$formula_nb)
  rhs_zi <- deparse_rhs(object$formulas$formula_zi)

if(zinb_comparison){
  f_zinb <- as.formula(paste(dv_f, '~', rhs_nb, '|', rhs_zi))
}
if(zip_comparison){
  f_zip <- as.formula(paste(dv_f, '~', rhs_nb, '|', rhs_zi))
}
  if(nb_comparison){
    if(!is.null(init_theta)){
  full_nb <- cmp_fit(MASS::glm.nb, list(formula = object$formulas$formula_nb, data = object$data$data, init.theta = init_theta))
    }else{
      full_nb <- cmp_fit(MASS::glm.nb, list(formula = object$formulas$formula_nb, data = object$data$data))
    }
  }
  if(zinb_comparison){
  full_zinb <- cmp_fit(pscl::zeroinfl, list(formula = f_zinb, data = object$data$data, dist = 'negbin'))
  }
  if(poisson_comparison){
  full_poisson <- cmp_fit(stats::glm, list(formula = object$formulas$formula_nb, data = object$data$data, family = stats::poisson()))
  }
  if(zip_comparison){
  full_zip <- cmp_fit(pscl::zeroinfl, list(formula = f_zip, data = object$data$data, dist = 'poisson'))
  }

  if(winsorize){
  #data_winsor <- object$data$data %>% dplyr::mutate(osvAll = dplyr::case_when(osvAll > sort(object$data$data$osvAll,decreasing=T)[cutoff_value] ~ sort(object$data$data$osvAll,decreasing=T)[cutoff_value],
  #                                                         T ~ osvAll))
    data_winsor <- object$data$data
    data_winsor[[dv_f]][which(data_winsor[dv_f]>sort(dplyr::pull(data_winsor[dv_f]),decreasing = T)[cutoff_value])] <- sort(dplyr::pull(data_winsor[dv_f]),decreasing = T)[cutoff_value]
  if(nb_comparison){
    if(!is.null(init_theta)){
  full_nb_winsor <- cmp_fit(MASS::glm.nb, list(formula = object$formulas$formula_nb, data = data_winsor, init.theta = init_theta))
    }else{
      full_nb_winsor <- cmp_fit(MASS::glm.nb, list(formula = object$formulas$formula_nb, data = data_winsor))
    }
  }
  if(zinb_comparison){
  full_zinb_winsor <- cmp_fit(pscl::zeroinfl, list(formula = f_zinb, data = data_winsor, dist = 'negbin'))
  }
  if(poisson_comparison){
  full_poisson_winsor <- cmp_fit(stats::glm, list(formula = object$formulas$formula_nb, data = data_winsor, family = stats::poisson()))
  }
  if(zip_comparison){
  full_zip_winsor <- cmp_fit(pscl::zeroinfl, list(formula = f_zip, data = data_winsor, dist = 'poisson'))
  }
  }
  if(razorize){
  # audit0.10 §1.5: keep_rows is reused below so boot_refit_one() can map a
  # bootstrap resample (indexed into the full data) onto data_razor's own
  # (fewer, renumbered) rows.
  y_full_for_razor <- dplyr::pull(object$data$data[dv_f])
  keep_rows <- which(y_full_for_razor < sort(y_full_for_razor,decreasing=T)[cutoff_value])
  data_razor <- object$data$data[keep_rows, ]
  if(nb_comparison){
    if(!is.null(init_theta)){
  full_nb_razor <- cmp_fit(MASS::glm.nb, list(formula = object$formulas$formula_nb, data = data_razor, init.theta = init_theta))
    }else{
      full_nb_razor <- cmp_fit(MASS::glm.nb, list(formula = object$formulas$formula_nb, data = data_razor))

    }
  }
  if(zinb_comparison){
  full_zinb_razor <- cmp_fit(pscl::zeroinfl, list(formula = f_zinb, data = data_razor, dist = 'negbin'))
  }
  if(poisson_comparison){
  full_poisson_razor <- cmp_fit(stats::glm, list(formula = object$formulas$formula_nb, data = data_razor, family = stats::poisson()))
  }
  if(zip_comparison){
  full_zip_razor <- cmp_fit(pscl::zeroinfl, list(formula = f_zip, data = data_razor, dist = 'poisson'))
  }
  }
  # One spec per family member: its (winsorised / razorised) data set, the
  # full-data fit computed above, and just enough to refit it on a resample
  # (no closures over `object`, so parallel workers get a small spec).
  fmls <- object$formulas
  # audit0.10 §1.5: the original (unwinsorised, unrazorised) response, so OOB
  # error is always computed against the real outcome regardless of what
  # `data` holds. For razorised specs it is subset to keep_rows below so it
  # stays aligned with data_razor's own row positions.
  y_orig_full <- dplyr::pull(object$data$data[dv_f])
  mk <- function(type, cls, data, full, keep = NULL) {
    list(class = cls, type = type, data = data, full = full,
         formulas = fmls, f_zinb = if (type == "zinb") f_zinb else NULL,
         f_zip = if (type == "zip") f_zip else NULL,
         has_init_theta = !is.null(init_theta), init_theta = init_theta,
         keep = keep, weights_col = wcol,
         y_orig = if (is.null(keep)) y_orig_full else y_orig_full[keep])
  }
  specs <- list()
  if (nb_comparison)      specs$nb      <- mk("nb",      "nbboot",      object$data$data, full_nb)
  if (zinb_comparison)    specs$zinb    <- mk("zinb",    "zinbboot",    object$data$data, full_zinb)
  if (poisson_comparison) specs$poisson <- mk("poisson", "poissonboot", object$data$data, full_poisson)
  if (zip_comparison)     specs$zip     <- mk("zip",     "zipboot",     object$data$data, full_zip)
  if (winsorize && nb_comparison)      specs$nb_winsor      <- mk("nb",      "nbboot",      data_winsor, full_nb_winsor)
  if (winsorize && zinb_comparison)    specs$zinb_winsor    <- mk("zinb",    "zinbboot",    data_winsor, full_zinb_winsor)
  if (winsorize && poisson_comparison) specs$poisson_winsor <- mk("poisson", "poissonboot", data_winsor, full_poisson_winsor)
  if (winsorize && zip_comparison)     specs$zip_winsor     <- mk("zip",     "zipboot",     data_winsor, full_zip_winsor)
  if (razorize && nb_comparison)      specs$nb_razor      <- mk("nb",      "nbboot",      data_razor, full_nb_razor, keep = keep_rows)
  if (razorize && zinb_comparison)    specs$zinb_razor    <- mk("zinb",    "zinbboot",    data_razor, full_zinb_razor, keep = keep_rows)
  if (razorize && poisson_comparison) specs$poisson_razor <- mk("poisson", "poissonboot", data_razor, full_poisson_razor, keep = keep_rows)
  if (razorize && zip_comparison)     specs$zip_razor     <- mk("zip",     "zipboot",     data_razor, full_zip_razor, keep = keep_rows)

  fam <- evinf_with_plan(multicore, ncores, {
    boot_refit_family(specs, object$bootstraps, object$boot_seeds[[1]])
  })

  out <- c(list(model = object), fam)
  # Backwards-compatible alias for the pre-0.9.4 name of the first slot.
  out$evzinb <- object
  class(out) <- 'evzinbcomp'
  return(out)
}

# Refit one bootstrap resample under one family spec (top-level so parallel
# workers do not receive the calling frame).
boot_refit_one <- function(sp, boot_id) {
  if (!is.null(sp$keep)) {
    # audit0.10 §1.5: boot_id indexes the full model data, but sp$data is
    # data_razor with fewer, renumbered rows. Keep only the drawn rows that
    # survived razorising and renumber them to data_razor's own positions;
    # inner_nb()/inner_zinb() then get the OOB set right "for free" via
    # ordinary negative indexing (data[-boot_id, ]).
    boot_id <- match(boot_id[boot_id %in% sp$keep], sp$keep)
  }
  b_stub <- list(boot_id = boot_id)
  if (sp$type == "nb") {
    if (sp$has_init_theta) {
      try(inner_nb(b_stub, sp$data, sp$formulas, sp$init_theta, sp$y_orig, sp$weights_col), silent = TRUE)
    } else {
      try(inner_nb(b_stub, sp$data, sp$formulas, y_orig = sp$y_orig, weights_col = sp$weights_col), silent = TRUE)
    }
  } else if (sp$type == "poisson") {
    try(inner_poisson(b_stub, sp$data, sp$formulas, sp$y_orig, sp$weights_col), silent = TRUE)
  } else if (sp$type == "zip") {
    try(inner_zip(b_stub, sp$data, sp$formulas, sp$f_zip, sp$y_orig, sp$weights_col), silent = TRUE)
  } else {
    try(inner_zinb(b_stub, sp$data, sp$formulas, sp$f_zinb, sp$y_orig, sp$weights_col), silent = TRUE)
  }
}

# Refit every family member (nb / zinb + winsorised / razorised variants) on
# every bootstrap resample, in one parallel pass over the (spec x bootstrap)
# grid (audit 4.1, round 5).
#
# @param specs Named list of family specs (see mk() in compare_models()).
# @param boots The \code{object$bootstraps} list (only \code{$boot_id} is used).
# @param seed Integer seed for evinf_pmap().
# @return A named list mirroring \code{specs}; each element is
#   \code{list(full_run, bootstraps)} with class \code{"nbboot"} / \code{"zinbboot"}.
# @noRd
boot_refit_family <- function(specs, boots, seed) {
  boot_ids <- purrr::map(boots, "boot_id")
  nb <- length(boot_ids)
  grid <- expand.grid(si = seq_along(specs), bi = seq_len(nb))

  flat <- evinf_pmap(
    seq_len(nrow(grid)),
    function(k, grid, specs, boot_ids)
      boot_refit_one(specs[[grid$si[k]]], boot_ids[[grid$bi[k]]]),
    grid = grid, specs = specs, boot_ids = boot_ids,
    seed = seed, label = "compare_models refit"
  )

  stats::setNames(lapply(seq_along(specs), function(si) {
    sp <- specs[[si]]
    fitted <- flat[grid$si == si]
    names(fitted) <- names(boots)
    fitted <- purrr::map(fitted, mr_inner)
    structure(list(full_run = sp$full, bootstraps = fitted), class = sp$class)
  }), names(specs))
}

# round8 0.7 (review §7): `data[-boot_id, ]` with `boot_id` of length zero
# returns zero rows (via `data[integer(0), ]`), not all of them -- a
# pathological but reachable case if razorising leaves no drawn row behind.
# Fail the same way a caught error from MASS::glm.nb()/pscl::zeroinfl()
# would, so every existing `inherits(b, "try-error")` check already handles it.
empty_boot_id_error <- function(fn_name) {
  try(stop(
    fn_name, "(): boot_id is empty for this resample (no rows survived ",
    "razorising), so data[-boot_id, ] would silently return zero rows ",
    "instead of all of them; skipping this replicate."
  ), silent = TRUE)
}

inner_nb <- function(bootstrap,data,formulas,init_theta,y_orig,weights_col=NULL){
  if (length(bootstrap$boot_id) == 0) {
    return(empty_boot_id_error("inner_nb"))
  }

  data_ib <- data[bootstrap$boot_id,]
  data_oob <- data[-bootstrap$boot_id,]
  # audit0.10 §1.5: the original (unwinsorised) outcome for these OOB rows,
  # not whatever `data` currently holds (winsorised values would make the
  # winsorised OOB error incomparable to evinf's).
  dv <- y_orig[-bootstrap$boot_id]
  # round10 0.2: weights_col names a column already present in data_ib
  # (resampling rows keeps every column, weights included); referenced by
  # name via do.call() rather than as a vector, see cmp_fit() above.
  args <- list(formula = formulas$formula_nb, data = data_ib)
  if (!missing(init_theta)) args$init.theta <- init_theta
  if (!is.null(weights_col)) args$weights <- as.name(weights_col)
  boot_nb <- try(do.call(MASS::glm.nb, args), silent = TRUE)
  if(!('try-error' %in% class(boot_nb))){
  boot_nb$oob_predictions <- exp(predict(boot_nb,newdata=data_oob))
  boot_nb$oob_rmse <- sqrt(mean((dv-boot_nb$oob_predictions)^2))
  boot_nb$oob_rmsle <- sqrt(mean((log1p(dv)-log1p(boot_nb$oob_predictions))^2))
  # audit 4.7: keep the per-bootstrap fit statistics for compare_fit().
  boot_nb$fit_stats <- c(logLik = as.numeric(stats::logLik(boot_nb)),
                         npar = boot_nb$rank + 1,
                         AIC = stats::AIC(boot_nb),
                         BIC = stats::BIC(boot_nb))
  boot_nb$model <- NULL
  boot_nb$y <- NULL
  boot_nb$linear.predictors <- NULL
  boot_nb$weights <- NULL
  boot_nb$residuals <- NULL
  boot_nb$fitted.values <- NULL
  boot_nb$effects <- NULL
  }
  return(boot_nb)
}

inner_zinb <- function(bootstrap,data,formulas,f_zinb,y_orig,weights_col=NULL){
  if (length(bootstrap$boot_id) == 0) {
    return(empty_boot_id_error("inner_zinb"))
  }

  data_ib <- data[bootstrap$boot_id,]
  data_oob <- data[-bootstrap$boot_id,]
  # audit0.10 §1.5: see inner_nb() -- original outcome, not `data`'s.
  dv <- y_orig[-bootstrap$boot_id]
  # Use the full two-part formula (count | zero), matching the full-sample fit
  # (audit 1.2); formulas$formula_zi alone would use the ZI regressors for both
  # parts.
  args <- list(formula = f_zinb, data = data_ib, dist = 'negbin')
  if (!is.null(weights_col)) args$weights <- as.name(weights_col)
  boot_zinb <- try(do.call(pscl::zeroinfl, args), silent = TRUE)
  if(!('try-error' %in% class(boot_zinb))){
    boot_zinb$oob_predictions <- predict(boot_zinb,newdata=data_oob)
    boot_zinb$oob_rmse <- sqrt(mean((dv-boot_zinb$oob_predictions)^2))
    boot_zinb$oob_rmsle <- sqrt(mean((log1p(dv)-log1p(boot_zinb$oob_predictions))^2))
    # audit 4.7: keep the per-bootstrap fit statistics for compare_fit().
    ll <- stats::logLik(boot_zinb)
    boot_zinb$fit_stats <- c(logLik = as.numeric(ll),
                             npar = attr(ll, "df"),
                             AIC = stats::AIC(boot_zinb),
                             BIC = stats::BIC(boot_zinb))
    boot_zinb$model <- NULL
    boot_zinb$y <- NULL
    boot_zinb$weights <- NULL
    boot_zinb$residuals <- NULL
    boot_zinb$fitted.values <- NULL
  }
  return(boot_zinb)
}

# round9 E.3 (audit §5.5): Poisson competitor baseline, mirroring inner_nb()
# above but via glm(family = poisson()) -- no dispersion parameter, so
# npar = rank (not rank + 1).
inner_poisson <- function(bootstrap,data,formulas,y_orig,weights_col=NULL){
  if (length(bootstrap$boot_id) == 0) {
    return(empty_boot_id_error("inner_poisson"))
  }

  data_ib <- data[bootstrap$boot_id,]
  data_oob <- data[-bootstrap$boot_id,]
  dv <- y_orig[-bootstrap$boot_id]
  args <- list(formula = formulas$formula_nb, data = data_ib, family = stats::poisson())
  if (!is.null(weights_col)) args$weights <- as.name(weights_col)
  boot_pois <- try(do.call(stats::glm, args), silent = TRUE)
  if(!('try-error' %in% class(boot_pois))){
  boot_pois$oob_predictions <- exp(predict(boot_pois,newdata=data_oob))
  boot_pois$oob_rmse <- sqrt(mean((dv-boot_pois$oob_predictions)^2))
  boot_pois$oob_rmsle <- sqrt(mean((log1p(dv)-log1p(boot_pois$oob_predictions))^2))
  boot_pois$fit_stats <- c(logLik = as.numeric(stats::logLik(boot_pois)),
                           npar = boot_pois$rank,
                           AIC = stats::AIC(boot_pois),
                           BIC = stats::BIC(boot_pois))
  boot_pois$model <- NULL
  boot_pois$y <- NULL
  boot_pois$linear.predictors <- NULL
  boot_pois$weights <- NULL
  boot_pois$residuals <- NULL
  boot_pois$fitted.values <- NULL
  boot_pois$effects <- NULL
  }
  return(boot_pois)
}

# round9 E.3: zero-inflated-Poisson competitor baseline, mirroring
# inner_zinb() above but via pscl::zeroinfl(dist = "poisson").
inner_zip <- function(bootstrap,data,formulas,f_zip,y_orig,weights_col=NULL){
  if (length(bootstrap$boot_id) == 0) {
    return(empty_boot_id_error("inner_zip"))
  }

  data_ib <- data[bootstrap$boot_id,]
  data_oob <- data[-bootstrap$boot_id,]
  dv <- y_orig[-bootstrap$boot_id]
  args <- list(formula = f_zip, data = data_ib, dist = 'poisson')
  if (!is.null(weights_col)) args$weights <- as.name(weights_col)
  boot_zip <- try(do.call(pscl::zeroinfl, args), silent = TRUE)
  if(!('try-error' %in% class(boot_zip))){
    boot_zip$oob_predictions <- predict(boot_zip,newdata=data_oob)
    boot_zip$oob_rmse <- sqrt(mean((dv-boot_zip$oob_predictions)^2))
    boot_zip$oob_rmsle <- sqrt(mean((log1p(dv)-log1p(boot_zip$oob_predictions))^2))
    ll <- stats::logLik(boot_zip)
    boot_zip$fit_stats <- c(logLik = as.numeric(ll),
                            npar = attr(ll, "df"),
                            AIC = stats::AIC(boot_zip),
                            BIC = stats::BIC(boot_zip))
    boot_zip$model <- NULL
    boot_zip$y <- NULL
    boot_zip$weights <- NULL
    boot_zip$residuals <- NULL
    boot_zip$fitted.values <- NULL
  }
  return(boot_zip)
}


model_remover <- function(obj){
  if('evzinb' %in% class(obj$full_run)){
    return(obj)
  }
  obj$bootstraps<-obj$bootstraps %>% purrr::map(mr_inner)

  return(obj)
}

mr_inner <- function(obj){
  if('try-error'%in%class(obj)){
    return(obj)
  }
  cl <- class(obj)
  # Keep 'qr' and 'terms': predict.nbboot()/predict.zinbboot() need them for
  # newdata predictions.
  obj <- obj[!(names(obj) %in% c('model', 'residuals', 'fitted.values',
                                 'weights', 'y', 'linear.predictors',
                                 'prior.weights'))]
  class(obj) <- cl
  return(obj)
}


# Function to obtain predicted probabilities for zeroinfl
prob_from_znb <- function(znb,
                          newdata=NULL){
if(is.null(newdata)){
  newdata <- znb$model[,-1]
}
    if(!('matrix' %in% class(newdata))){
      forms <- znb_formula_extractor(znb$formula)
      if(is(forms, 'formula')){
      xdata_zi <- model.matrix(znb_formula_extractor(znb$formula),newdata)
      }else{
        xdata_zi <- model.matrix(forms$zi,newdata)[,-1]
      }
    }else{
  xdata_zi <- cbind(1,as.matrix(newdata))
    }
    
  

  pr_zc <- exp(xdata_zi %*% znb$coefficients$zero)/(1+exp(xdata_zi %*% znb$coefficients$zero))
  pr_count <- 1-exp(xdata_zi %*% znb$coefficients$zero)/(1+exp(xdata_zi %*% znb$coefficients$zero))
  out <- tibble::tibble(pr_zc=as.numeric(pr_zc),
                pr_count = as.numeric(pr_count))
  return(out)
}

count_from_znb <- function(znb,
                           newdata=NULL){
  if(is.null(newdata)){
    newdata <- znb$model[,-1]
  }
  if(!('matrix' %in% class(newdata))){
    forms <- znb_formula_extractor(znb$formula)
    if(is(forms,'formula')){
      xdata_nb <- model.matrix(znb_formula_extractor(znb$formula),newdata)
    }else{
      xdata_nb <- model.matrix(forms$nb,newdata)[,-1]
    }
  }else{
    xdata_nb <- cbind(1,as.matrix(newdata))
  }
  
  count <- exp(xdata_nb %*% znb$coefficients$count)
  return(count)
}

znb_formula_extractor <- function(formula){
  tmp <- as.character(formula)[3]
  strs <- trimws(strsplit(tmp, '|', fixed = TRUE)[[1]])

  if(identical(strs[1],strs[2]) | length(strs) == 1){
    return(as.formula(paste0('~',strs[1])))
  }else{
    return(list(nb = as.formula(paste0('~',strs[1])),
                zi = as.formula(paste0('~',strs[2]))))
  }
}


quantiles_zinb <- function(quantile,mu,theta,p_zero){
  qnbinom(ifelse(quantile-p_zero>0,(quantile-p_zero)/(1-p_zero),0),mu=mu,size=theta)
}

quantiles_from_zinb <- function(quantile,znb,
                                newdata = NULL){


  prbs <- prob_from_znb(znb,
                        newdata = newdata)

  cnts <- count_from_znb(znb,
                         newdata = newdata)
  
  out <- quantiles_zinb(quantile,mu=cnts,theta=znb$theta,p_zero=prbs$pr_zc)

  return(out)
}



quantiles_from_nb <- function(quantile,nb,
                              newdata = NULL){

    if(is.null(newdata)){
      out <- qnbinom(quantile,mu=exp(predict(nb)),size=nb$theta)
    }else{
      out <- qnbinom(quantile,mu=exp(predict(nb,newdata=newdata)),size=nb$theta)
    }
    return(out)
}

# round9 E.3: Poisson / ZIP counterparts of quantiles_from_nb()/
# quantiles_zinb()/quantiles_from_zinb() above -- qpois() has no dispersion
# argument, so these are otherwise identical. prob_from_znb()/count_from_znb()
# are reused as-is for ZIP: neither touches the count distribution's
# dispersion, only pscl::zeroinfl()'s zero/count coefficients, which are
# structured the same way regardless of dist.
quantiles_from_poisson <- function(quantile,pois,
                                   newdata = NULL){
  if(is.null(newdata)){
    out <- qpois(quantile,lambda=exp(predict(pois)))
  }else{
    out <- qpois(quantile,lambda=exp(predict(pois,newdata=newdata)))
  }
  return(out)
}

quantiles_zip <- function(quantile,mu,p_zero){
  qpois(ifelse(quantile-p_zero>0,(quantile-p_zero)/(1-p_zero),0),lambda=mu)
}

quantiles_from_zip <- function(quantile,zip,
                               newdata = NULL){
  prbs <- prob_from_znb(zip,
                        newdata = newdata)
  cnts <- count_from_znb(zip,
                         newdata = newdata)
  out <- quantiles_zip(quantile,mu=cnts,p_zero=prbs$pr_zc)
  return(out)
}


#' Prediction for zinbboot
#'
#' @param object a fitted zinbboot object
#' @param newdata Data to make predictions on
#' @param type What prediction should be computed? One of \code{"predicted"},
#'   \code{"counts"}, \code{"zi"}, \code{"count_state"}, \code{"states"},
#'   \code{"all"} or \code{"quantile"}. (A zero-inflated negative binomial has no
#'   extreme-value state, so \code{"evinf"} is not accepted.)
#' @param pred Prediction type, 'original', 'bootstrap_median', or 'bootstrap_mean'
#' @param quantile Quantile for quantile prediction
#' @param confint Should confidence intervals be created?
#' @param conf_level Confidence level when predicting with CIs
#' @param ... Not used
#' 
#' @importFrom rlang :=
#'
#' @return Predictions from zinbboot
#' @export
predict.zinbboot <- function(object,newdata=NULL, type = c('predicted','counts','zi','count_state','states','all', 'quantile'), pred = c('original','bootstrap_median','bootstrap_mean'),quantile=NULL,confint=FALSE, conf_level=0.9,...){
  
  pred <- match.arg(pred, c('original','bootstrap_median','bootstrap_mean'))
  
  type <- match.arg(type,c('predicted','counts','zi','count_state','states','all', 'quantile'))
  
  if(type %in% c('states','all') & confint){
    stop('Confidence interval prediction only available for vector outputs')
  }
  
  if(pred %in% c('bootstrap_median','bootstrap_mean') | confint){
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))    
    nboots <- length(object$bootstraps)
    if(is.null(newdata)){
      newdata <- object$full_run$model
    }
    prbs_boot <- purrr::map(object$bootstraps, function(b)
      dplyr::bind_cols(prob_from_znb(b, newdata = newdata),
                       tibble::tibble(id = 1:nrow(newdata))))
    cnts_boot <- purrr::map(object$bootstraps, function(b)
      dplyr::bind_cols(tibble::tibble(count = count_from_znb(b, newdata = newdata)),
                       tibble::tibble(id = 1:nrow(newdata))))

    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
        q_boot <- purrr::map(object$bootstraps, function(b)
          tibble::tibble(q = quantiles_from_zinb(quantile, b, newdata = newdata),
                         id = 1:nrow(newdata)))
      }else{
        q_boot <- NULL
      }
    }else{
      q_boot <- NULL
    }

    prediction_boot <- purrr::map(object$bootstraps, function(b)
      tibble::tibble(pred = predict(b, type = 'r', newdata = newdata),
                     id = 1:nrow(newdata)))

  }
  
  
  if(pred == 'original'){
    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
        q <- quantiles_from_zinb(quantile,object$full_run,newdata = newdata)
      }else{
        q <- NULL
      }
    }
    ## Estimate component probabilities for all individuals
    prbs <- prob_from_znb(object$full_run,
                             newdata = newdata)
    ## Estimate mu_nb for all individuals
    cnts <- tibble::tibble(count = as.numeric(count_from_znb(object$full_run,
                               newdata = newdata)))

    predicted <- if(is.null(newdata)){
      predict(object$full_run, type = 'r')
    }else{
      predict(object$full_run, type = 'r', newdata = newdata)
    }

  }else if(pred=='bootstrap_median'){
    prbs <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-"id")
    cnts <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-"id")
    
    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(median) %>% dplyr::select(-"id") %>% dplyr::pull(.data$q)
    }else{
      q <- NULL
    }
    
    predicted <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% dplyr::summarize(predicted = median(.data$pred)) %>% dplyr::pull(.data$predicted)
    
  }else if(pred=='bootstrap_mean'){
    warning('Bootstrapped mean predictions are experimental and may yield infinite values')
    prbs <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-"id")
    cnts <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-"id")
  
    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(mean) %>% dplyr::select(-"id") %>% dplyr::pull(.data$q)
    }else{
      q <- NULL
    }
    predicted <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% dplyr::summarize(predicted = mean(.data$pred)) %>% dplyr::pull(.data$predicted)
  
  }
  
  
  if(confint){
    if(type %in% c('states','all')){
      stop("Confidence interval prediction only available for vector predictions (not 'states' or 'all')")
    }
    qs <- c((1-conf_level)/2,1-(1-conf_level)/2)
    
    if(type == 'predicted'){
      ci <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$pred,qs[1]),
                         ci_ub = quantile(.data$pred,qs[2])) %>%
        dplyr::select(-"id")
      return(dplyr::bind_cols(tibble::tibble(predicted=predicted),ci))
    }
   
    if(type == 'counts'){
      ci <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$count,qs[1]),
                         ci_ub = quantile(.data$count,qs[2])) %>% dplyr::select(-"id")
      return(dplyr::bind_cols(tibble::tibble(count=cnts$count),ci))
      
    }
    if(type == 'zi'){
      ci <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$pr_zc,qs[1]),
                         ci_ub = quantile(.data$pr_zc,qs[2])) %>% dplyr::select(-"id")
      return(dplyr::bind_cols(tibble::tibble(pr_zc=prbs$pr_zc),ci))
    }
  
    if(type == 'count_state'){
      ci <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$pr_count,qs[1]),
                         ci_ub = quantile(.data$pr_count,qs[2])) %>% dplyr::select(-"id")
      return(dplyr::bind_cols(tibble::tibble(pr_count=prbs$pr_count),ci))
    }
    if(type == 'quantile'){
      warning('Confidence interval prediction with Quantiles may yield unstable results')
      ci <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$q,qs[1]),
                         ci_ub = quantile(.data$q,qs[2])) %>% dplyr::select(-"id")
      q_name <- paste0('q',100*quantile)
      return(dplyr::bind_cols(tibble::tibble(!!q_name:=q),ci))
    }
    
  }else{
    if(type == "predicted"){
      return(predicted)
    }
    if(type == "counts"){
      return(cnts$count)
    }
    if(type == "zi"){
      return(prbs$pr_zc)
    }
    if(type == "count_state"){
      return(prbs$pr_count)
    }
    if(type == 'states'){
      return(prbs)
    }
    if(type == 'quantile'){
      return(q)
    }
    if(type == 'all'){
      q_name <- paste0('q',100*quantile)
      return(dplyr::bind_cols(tibble::tibble(predicted = predicted,
                                             !!q_name := q,
                                             prbs,cnts)))
    }
  }

}

#' Prediction for zipboot
#'
#' @param object a fitted zipboot object
#' @param newdata Data to make predictions on
#' @param type What prediction should be computed? One of \code{"predicted"},
#'   \code{"counts"}, \code{"zi"}, \code{"count_state"}, \code{"states"},
#'   \code{"all"} or \code{"quantile"}. (A zero-inflated Poisson has no
#'   extreme-value state, so \code{"evinf"} is not accepted.)
#' @param pred Prediction type, 'original', 'bootstrap_median', or 'bootstrap_mean'
#' @param quantile Quantile for quantile prediction
#' @param confint Should confidence intervals be created?
#' @param conf_level Confidence level when predicting with CIs
#' @param ... Not used
#'
#' @importFrom rlang :=
#'
#' @return Predictions from zipboot
#' @export
predict.zipboot <- function(object,newdata=NULL, type = c('predicted','counts','zi','count_state','states','all', 'quantile'), pred = c('original','bootstrap_median','bootstrap_mean'),quantile=NULL,confint=FALSE, conf_level=0.9,...){

  pred <- match.arg(pred, c('original','bootstrap_median','bootstrap_mean'))

  type <- match.arg(type,c('predicted','counts','zi','count_state','states','all', 'quantile'))

  if(type %in% c('states','all') & confint){
    stop('Confidence interval prediction only available for vector outputs')
  }

  if(pred %in% c('bootstrap_median','bootstrap_mean') | confint){
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    nboots <- length(object$bootstraps)
    if(is.null(newdata)){
      newdata <- object$full_run$model
    }
    prbs_boot <- purrr::map(object$bootstraps, function(b)
      dplyr::bind_cols(prob_from_znb(b, newdata = newdata),
                       tibble::tibble(id = 1:nrow(newdata))))
    cnts_boot <- purrr::map(object$bootstraps, function(b)
      dplyr::bind_cols(tibble::tibble(count = count_from_znb(b, newdata = newdata)),
                       tibble::tibble(id = 1:nrow(newdata))))

    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
        q_boot <- purrr::map(object$bootstraps, function(b)
          tibble::tibble(q = quantiles_from_zip(quantile, b, newdata = newdata),
                         id = 1:nrow(newdata)))
      }else{
        q_boot <- NULL
      }
    }else{
      q_boot <- NULL
    }

    prediction_boot <- purrr::map(object$bootstraps, function(b)
      tibble::tibble(pred = predict(b, type = 'r', newdata = newdata),
                     id = 1:nrow(newdata)))

  }


  if(pred == 'original'){
    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
        q <- quantiles_from_zip(quantile,object$full_run,newdata = newdata)
      }else{
        q <- NULL
      }
    }
    ## Estimate component probabilities for all individuals
    prbs <- prob_from_znb(object$full_run,
                             newdata = newdata)
    ## Estimate mu_nb for all individuals
    cnts <- tibble::tibble(count = as.numeric(count_from_znb(object$full_run,
                               newdata = newdata)))

    predicted <- if(is.null(newdata)){
      predict(object$full_run, type = 'r')
    }else{
      predict(object$full_run, type = 'r', newdata = newdata)
    }

  }else if(pred=='bootstrap_median'){
    prbs <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-"id")
    cnts <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-"id")

    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(median) %>% dplyr::select(-"id") %>% dplyr::pull(.data$q)
    }else{
      q <- NULL
    }

    predicted <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% dplyr::summarize(predicted = median(.data$pred)) %>% dplyr::pull(.data$predicted)

  }else if(pred=='bootstrap_mean'){
    warning('Bootstrapped mean predictions are experimental and may yield infinite values')
    prbs <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-"id")
    cnts <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-"id")

    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(mean) %>% dplyr::select(-"id") %>% dplyr::pull(.data$q)
    }else{
      q <- NULL
    }
    predicted <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% dplyr::summarize(predicted = mean(.data$pred)) %>% dplyr::pull(.data$predicted)

  }


  if(confint){
    if(type %in% c('states','all')){
      stop("Confidence interval prediction only available for vector predictions (not 'states' or 'all')")
    }
    qs <- c((1-conf_level)/2,1-(1-conf_level)/2)

    if(type == 'predicted'){
      ci <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize(ci_lb = quantile(.data$pred,qs[1]),
                         ci_ub = quantile(.data$pred,qs[2])) %>%
        dplyr::select(-"id")
      return(dplyr::bind_cols(tibble::tibble(predicted=predicted),ci))
    }

    if(type == 'counts'){
      ci <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize(ci_lb = quantile(.data$count,qs[1]),
                         ci_ub = quantile(.data$count,qs[2])) %>% dplyr::select(-"id")
      return(dplyr::bind_cols(tibble::tibble(count=cnts$count),ci))

    }
    if(type == 'zi'){
      ci <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize(ci_lb = quantile(.data$pr_zc,qs[1]),
                         ci_ub = quantile(.data$pr_zc,qs[2])) %>% dplyr::select(-"id")
      return(dplyr::bind_cols(tibble::tibble(pr_zc=prbs$pr_zc),ci))
    }

    if(type == 'count_state'){
      ci <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize(ci_lb = quantile(.data$pr_count,qs[1]),
                         ci_ub = quantile(.data$pr_count,qs[2])) %>% dplyr::select(-"id")
      return(dplyr::bind_cols(tibble::tibble(pr_count=prbs$pr_count),ci))
    }
    if(type == 'quantile'){
      warning('Confidence interval prediction with Quantiles may yield unstable results')
      ci <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize(ci_lb = quantile(.data$q,qs[1]),
                         ci_ub = quantile(.data$q,qs[2])) %>% dplyr::select(-"id")
      q_name <- paste0('q',100*quantile)
      return(dplyr::bind_cols(tibble::tibble(!!q_name:=q),ci))
    }

  }else{
    if(type == "predicted"){
      return(predicted)
    }
    if(type == "counts"){
      return(cnts$count)
    }
    if(type == "zi"){
      return(prbs$pr_zc)
    }
    if(type == "count_state"){
      return(prbs$pr_count)
    }
    if(type == 'states'){
      return(prbs)
    }
    if(type == 'quantile'){
      return(q)
    }
    if(type == 'all'){
      q_name <- paste0('q',100*quantile)
      return(dplyr::bind_cols(tibble::tibble(predicted = predicted,
                                             !!q_name := q,
                                             prbs,cnts)))
    }
  }

}

#' Prediction for nbboot
#'
#' @param object a fitted nbboot object
#' @param newdata Data to make predictions on
#' @param type What prediction should be computed?
#' @param pred Prediction type, 'original', 'bootstrap_median', or 'bootstrap_mean'
#' @param quantile Quantile for quantile prediction
#' @param confint Should confidence intervals be created?
#' @param conf_level Confidence level when predicting with CIs
#' @param ... Not used
#'
#' @importFrom rlang :=
#'
#' @return Predictions from nbboot
#' @export
predict.nbboot <- function(object,newdata=NULL, type = c('predicted','all', 'quantile'), pred = c('original','bootstrap_median','bootstrap_mean'),quantile=NULL,confint=FALSE, conf_level=0.9,...){
  
  pred <- match.arg(pred, c('original','bootstrap_median','bootstrap_mean'))
  
  type <- match.arg(type,c('predicted','all', 'quantile'))
  
  if(type %in% c('states','all') & confint){
    stop('Confidence interval prediction only available for vector outputs')
  }
  
  if(pred %in% c('bootstrap_median','bootstrap_mean') | confint){
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))    
    nboots <- length(object$bootstraps)
    if(is.null(newdata)){
      newdata <- object$full_run$model
    }
  
    
    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
        q_boot <- purrr::map(object$bootstraps, function(b)
          tibble::tibble(q = quantiles_from_nb(quantile, b, newdata = newdata),
                         id = 1:nrow(newdata)))
      }else{
        q_boot <- NULL
      }
    }else{
      q_boot <- NULL
    }

    prediction_boot <- purrr::map(object$bootstraps, function(b)
      tibble::tibble(pred = predict(b, type = 'r', newdata = newdata),
                     id = 1:nrow(newdata)))
    
  }
  
  
  if(pred == 'original'){
    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
        q <- quantiles_from_nb(quantile,object$full_run,newdata = newdata)
      }else{
        q <- NULL
      }
    }
  
    ## Estimate mu_nb for all individuals

    predicted <- if(is.null(newdata)){
      predict(object$full_run, type = 'r')
    }else{
      predict(object$full_run, type = 'r', newdata = newdata)
    }

  }else if(pred=='bootstrap_median'){

    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(median) %>% dplyr::select(-"id") %>% dplyr::pull(.data$q)
    }else{
      q <- NULL
    }
    
    predicted <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% dplyr::summarize(predicted = median(.data$pred)) %>% dplyr::pull(.data$predicted)
    
  }else if(pred=='bootstrap_mean'){
    warning('Bootstrapped mean predictions are experimental and may yield infinite values')
  
    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(mean) %>% dplyr::select(-"id") %>% dplyr::pull(.data$q)
    }else{
      q <- NULL
    }
    predicted <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% dplyr::summarize(predicted = mean(.data$pred)) %>% dplyr::pull(.data$predicted)
    
  }
  
  
  if(confint){
    if(type %in% c('states','all')){
      stop("Confidence interval prediction only available for vector predictions (not 'states' or 'all')")
    }
    qs <- c((1-conf_level)/2,1-(1-conf_level)/2)
    
    if(type == 'predicted'){
      ci <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$pred,qs[1]),
                         ci_ub = quantile(.data$pred,qs[2])) %>%
        dplyr::select(-"id")
      return(dplyr::bind_cols(tibble::tibble(predicted=predicted),ci))
    }
    if(type == 'quantile'){
      warning('Confidence interval prediction with Quantiles may yield unstable results')
      ci <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$q,qs[1]),
                         ci_ub = quantile(.data$q,qs[2])) %>% dplyr::select(-"id")
      q_name <- paste0('q',100*quantile)
      return(dplyr::bind_cols(tibble::tibble(!!q_name:=q),ci))
    }
    
  }else{
    if(type == "predicted"){
      return(predicted)
    }
    if(type == 'quantile'){
      return(q)
    }
    if(type == 'all'){
      q_name <- paste0('q',100*quantile)
      return(dplyr::bind_cols(tibble::tibble(predicted = predicted,
                                             !!q_name := q)))
    }
  }

}

#' Prediction for poissonboot
#'
#' @param object a fitted poissonboot object
#' @param newdata Data to make predictions on
#' @param type What prediction should be computed?
#' @param pred Prediction type, 'original', 'bootstrap_median', or 'bootstrap_mean'
#' @param quantile Quantile for quantile prediction
#' @param confint Should confidence intervals be created?
#' @param conf_level Confidence level when predicting with CIs
#' @param ... Not used
#'
#' @importFrom rlang :=
#'
#' @return Predictions from poissonboot
#' @export
predict.poissonboot <- function(object,newdata=NULL, type = c('predicted','all', 'quantile'), pred = c('original','bootstrap_median','bootstrap_mean'),quantile=NULL,confint=FALSE, conf_level=0.9,...){

  pred <- match.arg(pred, c('original','bootstrap_median','bootstrap_mean'))

  type <- match.arg(type,c('predicted','all', 'quantile'))

  if(type %in% c('states','all') & confint){
    stop('Confidence interval prediction only available for vector outputs')
  }

  if(pred %in% c('bootstrap_median','bootstrap_mean') | confint){
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))
    nboots <- length(object$bootstraps)
    if(is.null(newdata)){
      newdata <- object$full_run$model
    }


    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
        q_boot <- purrr::map(object$bootstraps, function(b)
          tibble::tibble(q = quantiles_from_poisson(quantile, b, newdata = newdata),
                         id = 1:nrow(newdata)))
      }else{
        q_boot <- NULL
      }
    }else{
      q_boot <- NULL
    }

    prediction_boot <- purrr::map(object$bootstraps, function(b)
      tibble::tibble(pred = predict(b, type = 'r', newdata = newdata),
                     id = 1:nrow(newdata)))

  }


  if(pred == 'original'){
    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
        q <- quantiles_from_poisson(quantile,object$full_run,newdata = newdata)
      }else{
        q <- NULL
      }
    }

    ## Estimate mu_nb for all individuals

    predicted <- if(is.null(newdata)){
      predict(object$full_run, type = 'r')
    }else{
      predict(object$full_run, type = 'r', newdata = newdata)
    }

  }else if(pred=='bootstrap_median'){

    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(median) %>% dplyr::select(-"id") %>% dplyr::pull(.data$q)
    }else{
      q <- NULL
    }

    predicted <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% dplyr::summarize(predicted = median(.data$pred)) %>% dplyr::pull(.data$predicted)

  }else if(pred=='bootstrap_mean'){
    warning('Bootstrapped mean predictions are experimental and may yield infinite values')

    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(mean) %>% dplyr::select(-"id") %>% dplyr::pull(.data$q)
    }else{
      q <- NULL
    }
    predicted <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% dplyr::summarize(predicted = mean(.data$pred)) %>% dplyr::pull(.data$predicted)

  }


  if(confint){
    if(type %in% c('states','all')){
      stop("Confidence interval prediction only available for vector predictions (not 'states' or 'all')")
    }
    qs <- c((1-conf_level)/2,1-(1-conf_level)/2)

    if(type == 'predicted'){
      ci <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize(ci_lb = quantile(.data$pred,qs[1]),
                         ci_ub = quantile(.data$pred,qs[2])) %>%
        dplyr::select(-"id")
      return(dplyr::bind_cols(tibble::tibble(predicted=predicted),ci))
    }
    if(type == 'quantile'){
      warning('Confidence interval prediction with Quantiles may yield unstable results')
      ci <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize(ci_lb = quantile(.data$q,qs[1]),
                         ci_ub = quantile(.data$q,qs[2])) %>% dplyr::select(-"id")
      q_name <- paste0('q',100*quantile)
      return(dplyr::bind_cols(tibble::tibble(!!q_name:=q),ci))
    }

  }else{
    if(type == "predicted"){
      return(predicted)
    }
    if(type == 'quantile'){
      return(q)
    }
    if(type == 'all'){
      q_name <- paste0('q',100*quantile)
      return(dplyr::bind_cols(tibble::tibble(predicted = predicted,
                                             !!q_name := q)))
    }
  }

}
