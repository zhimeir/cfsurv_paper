#' Confidence interval based on Cox model
#'
#' Construct conformal predictive interval based on cox model
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

cox_based <- function(x,c,alpha,
                      data_fit,
                      data_calib,
                      type,
                      dist,
                      weight_calib,
                      weight_new,
                      ftol=.1,tol=.1){
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
  
  if(type == "quantile"){
    ## Fit the survival model
    xnames <- paste0("X",1:p)
    fmla <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
    mdl <- survreg(fmla,data=data_fit,dist=dist)

    ## The fitted quantile for the calibration data
    res <- predict(mdl,
                  newdata = data_calib,
                  type="quantile",
                  p=alpha)
    quant <-  res  
    score <- pmin(c,quant)-data_calib$censored_T
  
    ## The fitted quantile for the new data
    newdata <- data.frame(x)
    colnames(newdata) <- xnames
    res <- predict(mdl,
                  newdata = newdata,
                  type="quantile",
                  p=alpha)
    new_quant <-  res

    ## Compute the calibration term
    calib_term <- sapply(X=weight_new,get_calibration,score=score,
                        weight_calib=weight_calib,alpha=alpha)
    ## obtain final confidence interval
    lower_bnd <- pmin(new_quant,c)-calib_term
    }

  if(type == "percentile"){
    ## Fit the model for S(y)=p(min(T,c)>=y|X)
    xnames <- paste0("X",1:p)
    data_fit <- data_fit[data_fit$C>=c,]
    data_fit$censored_T <- pmin(data_fit$censored_T,c)
   
    surv_data_fit <- data_fit
    surv_data_fit$censored_T <- -surv_data_fit$censored_T
    fmla <- with(surv_data_fit,as.formula(paste("censored_T ~ ", paste(xnames, collapse= "+"))))
    if(p==1){
      capture.output(bw <- npcdistbw(fmla),file='NULL')
    }else{
      capture.output(bw <- npcdistbw(fmla,ftol=ftol,tol=tol),file='NULL')
    }

    surv_data_calib <- data_calib
    surv_data_calib$censored_T <- -surv_data_calib$censored_T
    score<- npcdist(bws=bw,newdata = surv_data_calib)$condist
    
    ## Obtain the calibration term
    calib_term <- sapply(X=weight_new,get_calibration,score=score,
                         weight_calib=weight_calib,alpha=alpha)

    ## Obtain the final confidence interval
    lower_bnd <- rep(0,len_x)
    newdata <- data.frame(x)
    colnames(newdata) <- xnames
    for(i in 1:len_x){
      time_candidate <- seq(0,c+2,by=.1)
      score_candidate <- sapply(-time_candidate,get_survival_fun,x = newdata[i,],bw=bw,xnames=xnames) 
      ind <- min(which(score_candidate<=calib_term[i]))
      if(ind>1){
        lower_bnd[i] <- time_candidate[ind-1]
      }else{
        lower_bnd[i] <- 0
      }
      }
    }
 
  lower_bnd <- pmax(lower_bnd,0)
  lower_bnd <- pmin(lower_bnd,c)
  return(lower_bnd)
}


get_survival_fun <- function(x,t,bw,xnames){
  input_data <- data.frame(x)
  colnames(input_data) <- xnames
  input_data <- cbind(input_data,censored_T=t)
  val<- npcdist(bws=bw,newdata=input_data)$condist
  return(val)
}
