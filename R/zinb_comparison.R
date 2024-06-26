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
#' @param zinb_comparison Should comparions be made with the zinb model?
#' @param multicore Logical: should multiple cores be used
#' @param ncores Number of cores if multicore is used
#'
#' @return A list with the original model as the first object and compared models as the following objects
#' @export
#'
#' @examples 
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 10, multicore = TRUE, ncores = 2)
#' compare_models(model)
#' 
compare_models <- function(object, nb_comparison = TRUE, zinb_comparison = TRUE, winsorize = FALSE, razorize = FALSE, cutoff_value=10, init_theta=NULL, multicore = FALSE, ncores=NULL){

  dv_f <- all.vars(object$formulas$formula_nb)[1]
  iv_nb <- all.vars(object$formulas$formula_nb)[-1]
  iv_zi <- all.vars(object$formulas$formula_zi)[-1]
  if(is.null(iv_zi)){
    iv_zi <- '1'
  }
  
  i <- 'temp_iter'
if(zinb_comparison){
  f_zinb <- as.formula(paste(dv_f,'~',paste(iv_nb,collapse = '+'),'|',paste(iv_zi,collapse = '+')))
}
  if(nb_comparison){
    if(!is.null(init_theta)){
  full_nb <- try(MASS::glm.nb(object$formulas$formula_nb,data = object$data$data,init.theta=init_theta))
    }else{
      full_nb <- try(MASS::glm.nb(object$formulas$formula_nb,data = object$data$data))
    }
  }
  if(zinb_comparison){
  full_zinb <- try(pscl::zeroinfl(f_zinb,data = object$data$data,dist = 'negbin'))
  }
  
  if(winsorize){
  #data_winsor <- object$data$data %>% dplyr::mutate(osvAll = dplyr::case_when(osvAll > sort(object$data$data$osvAll,decreasing=T)[cutoff_value] ~ sort(object$data$data$osvAll,decreasing=T)[cutoff_value],
  #                                                         T ~ osvAll))
    data_winsor <- object$data$data
    data_winsor[[dv_f]][which(data_winsor[dv_f]>sort(dplyr::pull(data_winsor[dv_f]),decreasing = T)[cutoff_value])] <- sort(dplyr::pull(data_winsor[dv_f]),decreasing = T)[cutoff_value]
  if(nb_comparison){
    if(!is.null(init_theta)){
  full_nb_winsor <- try(MASS::glm.nb(object$formulas$formula_nb,data = data_winsor,init.theta=init_theta))
    }else{
      full_nb_winsor <- try(MASS::glm.nb(object$formulas$formula_nb,data = data_winsor))
    }
  }
  if(zinb_comparison){
  full_zinb_winsor <- try(pscl::zeroinfl(f_zinb,data = data_winsor,dist = 'negbin'))
  }
  }
  if(razorize){
  data_razor <- object$data$data %>% dplyr::filter(!!dplyr::sym(dv_f) < sort(dplyr::pull(object$data$data[dv_f]),decreasing=T)[cutoff_value])
  if(nb_comparison){
    if(!is.null(init_theta)){
  full_nb_razor <- try(MASS::glm.nb(object$formulas$formula_nb,data = data_razor,init.theta=init_theta))
    }else{
      full_nb_razor <- try(MASS::glm.nb(object$formulas$formula_nb,data = data_razor))
      
    }
  }
  if(zinb_comparison){
  full_zinb_razor <- try(pscl::zeroinfl(f_zinb,data = data_razor,dist = 'negbin'))
  }
  }
  if(multicore){
    if(is.null(ncores)){
      ncores <- parallel::detectCores()-1
    }
    if(.Platform$OS.type == 'windows'){
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
    }else{
      doParallel::registerDoParallel(cores = ncores)
    }
  }else{
    doParallel::stopImplicitCluster()
  }
if(nb_comparison){
  if(!is.null(init_theta)){
    if(.Platform$OS.type == 'windows'){
  bootstraps_nb <- foreach::foreach(i = 1:length(object$bootstraps),.packages = "evinf") %dopar%
    try(inner_nb(object$bootstraps[[i]],data = object$data$data,formulas=object$formulas,init_theta=init_theta))
    }else{
      bootstraps_nb <- foreach::foreach(i = 1:length(object$bootstraps)) %dopar%
        try(inner_nb(object$bootstraps[[i]],data = object$data$data,formulas=object$formulas,init_theta=init_theta))
    }
  }else{
    if(.Platform$OS.type == 'windows'){
    bootstraps_nb <- foreach::foreach(i = 1:length(object$bootstraps),.packages = "evinf") %dopar%
      try(inner_nb(object$bootstraps[[i]],data = object$data$data,formulas=object$formulas))
    }else{
      bootstraps_nb <- foreach::foreach(i = 1:length(object$bootstraps)) %dopar%
        try(inner_nb(object$bootstraps[[i]],data = object$data$data,formulas=object$formulas))
    }
  }
  names(bootstraps_nb) <- names(object$bootstraps)
  bootstraps_nb <- bootstraps_nb %>% purrr::map(mr_inner)
}
  if(zinb_comparison){
    if(.Platform$OS.type == 'windows'){
  bootstraps_zinb <- foreach::foreach(i = 1:length(object$bootstraps),.packages = "evinf") %dopar%
    try(inner_zinb(object$bootstraps[[i]],data = object$data$data,formulas = object$formulas))
    }else{
      bootstraps_zinb <- foreach::foreach(i = 1:length(object$bootstraps)) %dopar%
        try(inner_zinb(object$bootstraps[[i]],data = object$data$data,formulas = object$formulas))
    }
  names(bootstraps_zinb) <- names(object$bootstraps)
  bootstraps_zinb <- bootstraps_zinb %>% purrr::map(mr_inner)
  }
  if(winsorize){
  if(nb_comparison){
    if(!is.null(init_theta)){
      if(.Platform$OS.type == 'windows'){
    bootstraps_nb_winsor <- foreach::foreach(i = 1:length(object$bootstraps),.packages ="evinf") %dopar%
    try(inner_nb(object$bootstraps[[i]],data = data_winsor,formulas=object$formulas,init_theta=init_theta))
      }else{
        bootstraps_nb_winsor <- foreach::foreach(i = 1:length(object$bootstraps)) %dopar%
          try(inner_nb(object$bootstraps[[i]],data = data_winsor,formulas=object$formulas,init_theta=init_theta))
      }
    }else{
      if(.Platform$OS.type == 'windows'){
      bootstraps_nb_winsor <- foreach::foreach(i = 1:length(object$bootstraps),.packages = "evinf") %dopar%
        try(inner_nb(object$bootstraps[[i]],data = data_winsor,formulas=object$formulas))
      }else{
        bootstraps_nb_winsor <- foreach::foreach(i = 1:length(object$bootstraps)) %dopar%
          try(inner_nb(object$bootstraps[[i]],data = data_winsor,formulas=object$formulas))
      }
    }
  names(bootstraps_nb_winsor) <- names(object$bootstraps)
  bootstraps_nb_winsor<-bootstraps_nb_winsor %>% purrr::map(mr_inner)
  }
    if(zinb_comparison){
      if(.Platform$OS.type == 'windows'){
  bootstraps_zinb_winsor <- foreach::foreach(i = 1:length(object$bootstraps),.packages = "evinf") %dopar%
    try(inner_zinb(object$bootstraps[[i]],data = data_winsor,formulas = object$formulas))
      }else{
        bootstraps_zinb_winsor <- foreach::foreach(i = 1:length(object$bootstraps)) %dopar%
          try(inner_zinb(object$bootstraps[[i]],data = data_winsor,formulas = object$formulas))
      }
  names(bootstraps_zinb_winsor) <- names(object$bootstraps)
  bootstraps_zinb_winsor <- bootstraps_zinb_winsor %>% purrr::map(mr_inner)
    }
  }
if(razorize){
  if(nb_comparison){
    if(!is.null(init_theta)){
      if(.Platform$OS.type == 'windows'){
  bootstraps_nb_razor <- foreach::foreach(i = 1:length(object$bootstraps),.package = "evinf") %dopar%
    try(inner_nb(object$bootstraps[[i]],data = data_razor,formulas = object$formulas,init_theta=init_theta))
      }else{
        bootstraps_nb_razor <- foreach::foreach(i = 1:length(object$bootstraps)) %dopar%
          try(inner_nb(object$bootstraps[[i]],data = data_razor,formulas = object$formulas,init_theta=init_theta))
      }
    }else{
      if(.Platform$OS.type == 'windows'){
      bootstraps_nb_razor <- foreach::foreach(i = 1:length(object$bootstraps),.package = "evinf") %dopar%
        try(inner_nb(object$bootstraps[[i]],data = data_razor,formulas = object$formulas))
      }else{
        bootstraps_nb_razor <- foreach::foreach(i = 1:length(object$bootstraps)) %dopar%
          try(inner_nb(object$bootstraps[[i]],data = data_razor,formulas = object$formulas))
      }
    }
  names(bootstraps_nb_razor) <- names(object$bootstraps)
  bootstraps_nb_razor<- bootstraps_nb_razor %>% purrr::map(mr_inner)
  }
  if(zinb_comparison){
    if(.Platform$OS.type == 'windows'){
  bootstraps_zinb_razor <- foreach::foreach(i = 1:length(object$bootstraps),.package = "evinf") %dopar%
    try(inner_zinb(object$bootstraps[[i]],data = data_razor,formulas = object$formulas))
    }else{
      bootstraps_zinb_razor <- foreach::foreach(i = 1:length(object$bootstraps)) %dopar%
        try(inner_zinb(object$bootstraps[[i]],data = data_razor,formulas = object$formulas))
    }
  names(bootstraps_zinb_razor) <- names(object$bootstraps)

  bootstraps_zinb_razor<- bootstraps_zinb_razor %>% purrr::map(mr_inner)
  }
}
  if(multicore){
    if(.Platform$OS.type == 'windows'){
    parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }else{
    doParallel::stopImplicitCluster()
    }
  }
if(nb_comparison){
  nb <- list(full_run = full_nb,
             bootstraps = bootstraps_nb)
  class(nb) <- 'nbboot'
}
  if(zinb_comparison){
  zinb <- list(full_run = full_zinb,
               bootstraps = bootstraps_zinb)
  class(zinb) <- 'zinbboot'
  }
if(winsorize){
  if(nb_comparison){
  nb_winsor <- list(full_run = full_nb_winsor,
             bootstraps = bootstraps_nb_winsor)
  class(nb_winsor) <- 'nbboot'
  }
  if(zinb_comparison){
  zinb_winsor <- list(full_run = full_zinb_winsor,
                      bootstraps = bootstraps_zinb_winsor)
  class(zinb_winsor) <- 'zinbboot'
  }
}
  
  if(razorize){
    if(nb_comparison){
  nb_razor <- list(full_run = full_nb_razor,
             bootstraps = bootstraps_nb_razor)


   class(nb_razor) <-'nbboot'
    }
    if(zinb_comparison){
  zinb_razor <- list(full_run = full_zinb_razor,
                   bootstraps = bootstraps_zinb_razor)
class(zinb_razor) <- 'zinbboot'
}
}
  
  out <- list()
  out$evzinb <- object
  if(nb_comparison){
  out$nb <- nb
  }
  if(zinb_comparison){
  out$zinb <- zinb
  }
  if(winsorize){
    if(nb_comparison){
  out$nb_winsor <- nb_winsor
    }
    if(zinb_comparison){
  out$zinb_winsor <- zinb_winsor
    }
  }
  if(razorize){
    if(nb_comparison){
  out$nb_razor <- nb_razor
    }
    if(zinb_comparison){
  out$zinb_razor <- zinb_razor
    }
  }
  class(out) <- 'evzinbcomp'
  return(out)
}

inner_nb <- function(bootstrap,data,formulas,init_theta){


  data_ib <- data[bootstrap$boot_id,]
  data_oob <- data[-bootstrap$boot_id,]
  dv <- model.response(model.frame(formulas$formula_nb,data_oob))
  boot_nb <- try(MASS::glm.nb(formulas$formula_nb,data = data_ib,init.theta=init_theta))
  if(!('try-error' %in% class(boot_nb))){
  boot_nb$oob_predictions <- exp(predict(boot_nb,newdata=data_oob))
  boot_nb$oob_rmse <- sqrt(mean((dv-boot_nb$oob_predictions)^2))
  boot_nb$oob_rmsle <- sqrt(mean((log1p(dv)-log1p(boot_nb$oob_predictions))^2))
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

inner_zinb <- function(bootstrap,data,formulas){
  data_ib <- data[bootstrap$boot_id,]
  data_oob <- data[-bootstrap$boot_id,]
  dv <- model.response(model.frame(formulas$formula_nb,data_oob))
  boot_zinb <- try(pscl::zeroinfl(formulas$formula_zi,data = data_ib,dist = 'negbin'))
  if(!('try-error' %in% class(boot_zinb))){
    boot_zinb$oob_predictions <- predict(boot_zinb,newdata=data_oob)
    boot_zinb$oob_rmse <- sqrt(mean((dv-boot_zinb$oob_predictions)^2))
    boot_zinb$oob_rmsle <- sqrt(mean((log1p(dv)-log1p(boot_zinb$oob_predictions))^2))
    boot_zinb$model <- NULL
    boot_zinb$y <- NULL
    boot_zinb$weights <- NULL
    boot_zinb$residuals <- NULL
    boot_zinb$fitted.values <- NULL
  }
  return(boot_zinb)
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
  obj[!(names(obj)%in%c('model',"residuals","fitted.values",'weights','y','linear.predictors','prior.weights','qr'))]
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
  strs <- stringr::str_split(tmp, ' \\| ') %>% purrr::reduce(c)

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


#' Prediction for zinbboot
#'
#' @param object a fitted zinbboot object
#' @param newdata Data to make predictions on
#' @param type What prediction should be computed?
#' @param pred Prediction type, 'original', 'bootstra_median', or 'bootstrap_mean'
#' @param quantile Quantile for quantile prediction
#' @param confint Should confidence intervals be created?
#' @param conf_level Confidence level when predicting with CIs
#' @param ... Not used
#' 
#' @importFrom rlang :=
#'
#' @return Predictions from zinbboot
predict.zinbboot <- function(object,newdata=NULL, type = c('predicted','counts','zi','evinf','count_state','states','all', 'quantile'), pred = c('original','bootstrap_median','bootstrap_mean'),quantile=NULL,confint=FALSE, conf_level=0.9,...){
  
  pred <- match.arg(pred, c('original','bootstrap_median','bootstrap_mean'))
  
  type <- match.arg(type,c('predicted','counts','zi','evinf','count_state','states','all', 'quantile'))
  i <- 'temp_iter'
  
  if(type %in% c('states','all') & confint){
    stop('Confidence interval prediction only available for vector outputs')
  }
  
  if(pred %in% c('bootstrap_median','bootstrap_mean') | confint){
    object$bootstraps <- object$bootstraps %>% purrr::discard(~'try-error' %in% class(.x))    
    nboots <- length(object$bootstraps)
    if(is.null(newdata)){
      newdata <- object$full_run$model
    }
    prbs_boot <- foreach::foreach(i = 1:nboots) %do%
      dplyr::bind_cols(prob_from_znb(object$bootstraps[[i]],
                                        newdata = newdata),tibble::tibble(id=1:nrow(newdata)))
    cnts_boot <- foreach::foreach(i = 1:nboots) %do%
      dplyr::bind_cols(tibble::tibble(count = count_from_znb(object$bootstraps[[i]],
                                          newdata = newdata)),tibble::tibble(id=1:nrow(newdata)))
    
    if(type %in% c('quantile','all')){
      if(type == 'quantile' & is.null(quantile)){
        stop('quantile must be provided for quantile prediction')
      }else if(!is.null(quantile)){
        q_boot <- foreach::foreach(i=1:nboots) %do%
          tibble::tibble(q=quantiles_from_zinb(quantile,object$bootstraps[[i]],newdata = newdata), id=1:nrow(newdata))
      }else{
        q_boot <- NULL
      }
    }else{
      q_boot <- NULL
    }
    
    prediction_boot <-  foreach::foreach(i = 1:nboots) %do%
      tibble::tibble(pred = predict(object$bootstraps[[i]],type = 'r', newdata=newdata),
                                    id = 1:nrow(newdata))
    
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
    cnts <- count_from_znb(object$full_run,
                               newdata = newdata)

    predicted <- predict(object$bootstraps[[i]],type = 'r', newdata=newdata)
    
  }else if(pred=='bootstrap_median'){
    prbs <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-.data$id)
    cnts <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(median) %>% dplyr::select(-.data$id)
    
    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(median) %>% dplyr::select(-.data$id) %>% dplyr::pull(.data$q)
    }else{
      q <- NULL
    }
    
    predicted <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% dplyr::summarize(predicted = median(.data$pred)) %>% dplyr::pull(.data$predicted)
    
  }else if(pred=='bootstrap_mean'){
    warning('Bootstrapped mean predictions are experimental and may yield infinite values')
    prbs <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-.data$id)
    cnts <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
      dplyr::summarize_all(mean) %>% dplyr::select(-.data$id)
  
    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(mean) %>% dplyr::select(-.data$id) %>% dplyr::pull(.data$q)
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
        dplyr::select(-.data$id)
      return(dplyr::bind_cols(tibble::tibble(predicted=predicted),ci))
    }
   
    if(type == 'counts'){
      ci <- cnts_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$count,qs[1]),
                         ci_ub = quantile(.data$count,qs[2])) %>% dplyr::select(-.data$id)
      return(dplyr::bind_cols(tibble::tibble(count=cnts$count),ci))
      
    }
    if(type == 'zi'){
      ci <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$pr_zc,qs[1]),
                         ci_ub = quantile(.data$pr_zc,qs[2])) %>% dplyr::select(-.data$id)
      return(dplyr::bind_cols(tibble::tibble(pr_zc=prbs$pr_zc),ci))
    }
  
    if(type == 'count_state'){
      ci <- prbs_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$pr_count,qs[1]),
                         ci_ub = quantile(.data$pr_count,qs[2])) %>% dplyr::select(-.data$id)
      return(dplyr::bind_cols(tibble::tibble(pr_count=prbs$pr_count),ci))
    }
    if(type == 'quantile'){
      warning('Confidence interval prediction with Quantiles may yield unstable results')
      ci <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$q,qs[1]),
                         ci_ub = quantile(.data$q,qs[2])) %>% dplyr::select(-.data$id)
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

predict.nbboot <- function(object,newdata=NULL, type = c('predicted','all', 'quantile'), pred = c('original','bootstrap_median','bootstrap_mean'),quantile=NULL,confint=F, conf_level=0.9,...){
  
  pred <- match.arg(pred, c('original','bootstrap_median','bootstrap_mean'))
  
  type <- match.arg(type,c('predicted','all', 'quantile'))
  i <- 'temp_iter'
  
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
        q_boot <- foreach::foreach(i=1:nboots) %do%
          tibble::tibble(q=quantiles_from_nb(quantile,object$bootstraps[[i]],newdata = newdata), id=1:nrow(newdata))
      }else{
        q_boot <- NULL
      }
    }else{
      q_boot <- NULL
    }
    
    prediction_boot <-  foreach::foreach(i = 1:nboots) %do%
      tibble::tibble(pred = predict(object$bootstraps[[i]],type = 'r', newdata=newdata),
                     id = 1:nrow(newdata))
    
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
    
    predicted <- predict(object$bootstraps[[i]],type = 'r', newdata=newdata)
    
  }else if(pred=='bootstrap_median'){
 
    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(median) %>% dplyr::select(-.data$id) %>% dplyr::pull(.data$q)
    }else{
      q <- NULL
    }
    
    predicted <- prediction_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% dplyr::summarize(predicted = median(.data$pred)) %>% dplyr::pull(.data$predicted)
    
  }else if(pred=='bootstrap_mean'){
    warning('Bootstrapped mean predictions are experimental and may yield infinite values')
  
    if(!is.null(q_boot)){
      q <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>%
        dplyr::summarize_all(mean) %>% dplyr::select(-.data$id) %>% dplyr::pull(.data$q)
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
        dplyr::select(-.data$id)
      return(dplyr::bind_cols(tibble::tibble(predicted=predicted),ci))
    }
    if(type == 'quantile'){
      warning('Confidence interval prediction with Quantiles may yield unstable results')
      ci <- q_boot %>% dplyr::bind_rows() %>% dplyr::group_by(.data$id) %>% 
        dplyr::summarize(ci_lb = quantile(.data$q,qs[1]),
                         ci_ub = quantile(.data$q,qs[2])) %>% dplyr::select(-.data$id)
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
