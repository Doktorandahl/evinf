
#' Out of bag predictive performance of EVZINB and EVINB models
#'
#' @param object A fitted evzinb or evinb with bootstraps on which to conduct out-of-bag evaluation
#' @param predict_type What type of prediction should be made? Harmonic mean, or exp(log(prediction))?
#' @param metric What metric should be used for the out of bag evaluation? Default options include rmsle, rmse, mse, and mae. Can also take a user supplied function of the form function(y_pred,y_true) which returns a single value
#'
#' @return A vector of oob evaluation metrics of the length of the number of bootstraps in the evzinb/evinb object.
#' @export
#'
#' @examples 
#' data(genevzinb2)
#' model <- evzinb(y~x1+x2+x3,data=genevzinb2, n_bootstraps = 10)
#' oob_evaluation(model)
#' 
oob_evaluation <- function(object,predict_type = c('harmonic','explog'),
                           metric = c('rmsle','rmse','mse','mae')){
  i <- 'temp_iter'
  
  if(is.character(metric)){
    metric <- match.arg(metric, c('rmsle','rmse','mse','mae'))
    ev_metric <- switch(metric,
                        rmsle = rmsle_metric,
                        rmse = rmse_metric,
                        mse = mse_metric,
                        mae = mae_metric)
  }else{
    ev_metric <- metric
  }
  predict_type <- match.arg(predict_type, c('harmonic','explog'))
  
  evals <- foreach::foreach(i = 1:length(object$bootstraps)) %do%
    try(oob_inner(object$bootstraps[[i]],object$data,predict_type,ev_metric,model_type = class(object)))
  
  evals <- purrr::map(evals,err2na)
  evals <- purrr::reduce(evals,c)
  
  return(evals)
}


oob_inner <- function(bootstrap,data,predict_type,ev_metric,model_type = c('evzinb','evinb')){
  
  oob_data <- data$data[-bootstrap$boot_id,]
if(model_type=="evzinb"){
  predictions <- predict.evzinb(bootstrap,newdata=oob_data,type=predict_type)
}else{
  predictions <- predict.evinb(bootstrap,newdata=oob_data,type=predict_type)
}
  ev_metric(predictions,data$y[-bootstrap$boot_id])
  
}


err2na <- function(x){
  if(is.null(x) | is(x,'try-error')){
    return(NA)
  }else{
    return(x)
  }
}

# Base-R replacements for the handful of MLmetrics one-liners the package used.
# Signature matches MLmetrics: function(y_pred, y_true).
rmsle_metric <- function(y_pred, y_true) sqrt(mean((log1p(y_pred) - log1p(y_true))^2))
rmse_metric <- function(y_pred, y_true) sqrt(mean((y_pred - y_true)^2))
mse_metric <- function(y_pred, y_true) mean((y_pred - y_true)^2)
mae_metric <- function(y_pred, y_true) mean(abs(y_pred - y_true))


  
