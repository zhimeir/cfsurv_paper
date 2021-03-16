#' Confidence interval based on Powell (1986)
#'
#' Construct conformal predictive interval based on cox model
#'
#' @param x a vector of the covariate of the test data.
#' @param r the censoring time of the test data.
#' @param alpha a number betweeo 0 and 1, specifying the miscaverage rate.
#' @param data_fit a data frame, containing the training data.
#' @param data_calib a data frame, containing the calibration data.
#' @param type either "marginal" or "local". Determines the type of confidence interval.
#'
#' @return low_ci a value of the lower bound for the survival time of the test point.
#' @return includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
#'
#' @family model
#'
#' @export

pow_based <- function(x,c,alpha,
                      data_fit,
                      data_calib,
                      weight_calib,
                      weight_new){
  ## Check the dimensionality
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }

  ## Keep only the data points with C>=c and transform min(T,C) to min(T,c) 
  weight_calib <- weight_calib[data_calib$C>=c]
  data_calib <- data_calib[data_calib$C>=c,]
  data_calib$censored_T <- pmin(data_calib$censored_T,c)

  ## Fit the quantile regression model (Powell)
  xnames <- paste0("X",1:p)
  fmla <- as.formula(paste("Curv(censored_T,C,ctype=\"right\") ~ ", paste(xnames, collapse= "+")))
  mdl <- crq(fmla,data=data_fit,method = "Pow",taus = alpha)
  res <- predict(mdl,
        newdata = data_calib,
        type = "quantile")
  quant <- res

  newdata <- data.frame(x)
  colnames(newdata) <- xnames
  res <- predict(mdl,
        newdata = newdata,
        type = "quantile")
  new_quant <- res

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
