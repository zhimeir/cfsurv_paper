#' Confidence interval based on random forest
#'
#' Construct a conformal predictive confidence interval based on random forest
#'
#' @param x a vector of the covariate of the test data.
#' @param r the censoring time of the test data.
#' @param alpha a number betweeo 0 and 1, specifying the miscaverage rate.
#' @param data_fit a data frame, containing the training data.
#' @param data_calib a data frame, containing the calibration data.
#' @param type either "marginal" or "local". Determines the type of confidence interval.
#' @param dist The distribution of T used in the cox model.
#'
#' @return low_ci a value of the lower bound for the survival time of the test point.
#' @return includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
#'
#' @family model
#'
#' @export

rf_based <- function(x,c,alpha,
                     data_fit,
                     data_calib,
                     weight_calib,
                     weight_new){
  
  ## Keep only the data points with C>=c and transform min(T,C) to min(T,c) 
  weight_calib <- weight_calib[data_calib$C>=c]
  data_calib <- data_calib[data_calib$C>=c,]
  data_calib$censored_T <- pmin(data_calib$censored_T,c)
  
  ## Parameters
  n_calib <- dim(data_calib)[1]
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
    xnames <- paste0("X",1:p)
    data_test <- c(data_calib$X1,x)
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
    xnames <- paste0("X",1:p)
    data_test <- rbind(data_calib[,colnames(data_calib)%in%xnames],x)
  }


  ## Fit the model
  ntree <- 1000
  nodesize <- 80

  data_test <- data.frame(data_test)
  names(data_test) <- xnames
  fmla <- as.formula(paste("censored_T ~ ",paste(xnames,collapse="+")))
  mdl <- crf.km(fmla, ntree = ntree, 
                 nodesize = nodesize,
                 data_train = data_fit[,names(data_fit)%in%
                                       c(xnames,"censored_T","event")], 
                 data_test = data_test, 
                 yname = 'censored_T', 
                 iname = 'event',
                 tau = alpha,
                 method = "grf")
  quant <- mdl$predicted[1:n_calib]
  new_quant <- tail(mdl$predicted,-n_calib)
  score <- pmin(c,quant)-data_calib$censored_T
  
  ## Compute the calibration term
  calib_term <- sapply(X=weight_new,get_calibration,score=score,
                      weight_calib=weight_calib,alpha=alpha)
  ## obtain final confidence interval
  lower_bnd <- pmin(new_quant,c)-calib_term
  lower_bnd <- pmax(lower_bnd,0)
  lower_bnd <- pmin(lower_bnd,c)
  return(lower_bnd)
}
