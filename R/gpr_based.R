#' Confidence interval based on Gaussian process regression
#'
#' Construct conformal predictive interval based on distributional boosting
#'
#' @param x a vector of the covariate of the test data.
#' @param c the censoring time of the test data.
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

gpr_based <- function(x,c,alpha,
                      data_fit,
                      data_calib,
                      weight_calib,
                      weight_new){
  ## Check the dimensionality of the input
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }
  
  ## Keep only the data points with C>=c
  ## Transform min(T,C) to min(T,c) 
  weight_calib <- weight_calib[data_calib$C>=c]
  data_calib <- data_calib[data_calib$C>=c,]
  data_calib$censored_T <- pmin(data_calib$censored_T,c)
  
  ## Fit the model for S(y)=p(min(T,c)>=y|X)
  xnames <- paste0("X",1:p)
  data_fit <- data_fit[data_fit$C>=c,]
  data_fit$censored_T <- pmin(data_fit$censored_T,c)

  ## Moving on, use surv_data_fit  
  surv_data_fit <- data_fit
  surv_data_fit$censored_T <- -surv_data_fit$censored_T
  gpr_mdl <- GauPro(X = as.matrix(surv_data_fit[,names(surv_data_fit) %in% xnames]),
                    Z = surv_data_fit$censored_T, D = p,
                    type = "Gauss")

  ## obtain the score of the calibration data
  surv_data_calib <- data_calib
  surv_data_calib$censored_T <- -surv_data_calib$censored_T
  mean_calib <- gpr_mdl$predict(surv_data_calib[,
                                names(surv_data_calib) %in% xnames])
  sd_calib <- gpr_mdl$predict(surv_data_calib[,
                              names(surv_data_calib) %in% xnames],
                              se.fit = TRUE)$se

  score <- pnorm((surv_data_calib$censored_T - mean_calib) / sd_calib)
 
  ## Obtain the calibration term
  calib_term <- sapply(X=weight_new,get_calibration,score=score,
                         weight_calib=weight_calib,alpha=alpha)

  ## Obtain the final confidence interval
  lower_bnd <- rep(0,len_x)
  mean_new <- gpr_mdl$predict(x)
  sd_new <- gpr_mdl$predict(x, se.fit = TRUE)$se
  for(i in 1:len_x){
    lower_bnd[i] <- -(mean_new[i] + sd_new[i] * qnorm(calib_term[i]))
  }
 
  lower_bnd <- pmax(lower_bnd,0)
  lower_bnd <- pmin(lower_bnd,c)
  return(lower_bnd)
}



