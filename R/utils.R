#' Fitting the censoring probability P(C>=c|X)
#'
#' @export

censoring_prob <- function(fit, calib, test=NULL,
                           method = "distBoost",
                           xnames, c,
                           ftol=.1, tol=.1, n.tree = 40){

  p <- length(xnames)
  if(method == "np"){
    ## Fitting P(-C<=-c_0|X) (since P(C>=c_0|X)=P(-C<=-c_0|X))
    fit$C <- -fit$C
    fmla <- with(fit,as.formula(paste("C ~ ", paste(xnames, collapse= "+"))))
    if(length(xnames)==1){
      capture.output(bw <- npcdistbw(fmla),file =NULL)
    }else{
      capture.output(bw <- npcdistbw(fmla,ftol=ftol,tol=tol),file =NULL)
    }

    ## Computing censoring scores for the calibration data
    newdata_calib <- calib
    newdata_calib$C <- -c
    pr_calib<- npcdist(bws=bw,newdata = newdata_calib)$condist

    ## Computing the censoring scores for the test data
    if(!is.null(test)){
      newdata <- cbind(test,C=-c)
      newdata <- data.frame(newdata)
      colnames(newdata) <- c(xnames,"C")
      pr_new <- npcdist(bws=bw,newdata=newdata)$condist
    }else{pr_new=NULL}
  }

  if(method == "distBoost"){
    ## Fitting P(-C<=-c_0|X) (since P(C>=c_0|X)=P(-C<=-c_0|X))
    fit$C <- -fit$C
    fmla <- with(fit,as.formula(paste("C ~ ", paste(xnames, collapse= "+"))))
    gbm_mdl <- gbm(fmla,data=fit,distribution="gaussian", n.tree = n.tree)
    median_fit<- predict(object=gbm_mdl,newdata = fit)
    res_fit <- fit$C-median_fit
    resamp_fit <- median_fit + res_fit[sample.int(dim(fit)[1])]
    if(p==1){
      xdf <- data.frame(X1=fit[,colnames(fit)%in%xnames],
                    X2=rep(1,dim(fit)[1]))
     }else{
      xdf <- fit[,colnames(fit)%in%xnames]
     }
    mdlrb <- modtrast(xdf,fit$C,resamp_fit,min.node=200)
    
    ## Computing the censoring scores for the fitting data
    pr_fit <- rep(NA, dim(fit)[1])
    for(i in 1:length(pr_fit)){
        pr_fit[i] <- distBoost_cdf(mdlrb,xdf[i,],median_fit[i],-c,res_fit)
    }

    ## Computing the censoring scores for the calibration data
    pr_calib <- rep(NA, dim(calib)[1])
    median_calib<- predict(object=gbm_mdl,newdata = calib)
    if(p==1){
      xdf <- data.frame(X1=calib[,colnames(calib)%in%xnames],
                    X2=rep(1,dim(calib)[1]))
     }else{
      xdf <- calib[,colnames(calib)%in%xnames]
     }

    for(i in 1:length(pr_calib)){
      pr_calib[i] <- distBoost_cdf(mdlrb,xdf[i,],median_calib[i],
                               -c,res_fit)
    }

    ## Computing the censoring scores for the test data
    if(!is.null(test)){
      newdata <- data.frame(test)
      colnames(newdata) <- xnames
      median_test<- predict(object=gbm_mdl,newdata = newdata)
      if(p==1){
        xdf <- data.frame(X1=test,
                    X2=rep(1,length(test)))
        n_new <- length(test)
      }else{
        xdf <- test
        n_new <- dim(test)[1]
      }
      pr_new <- rep(NA, n_new)
      for(i in 1:n_new){
          pr_new[i] <- distBoost_cdf(mdlrb,xdf[i,],median_test[i],-c,res_fit)
      }
    }else{pr_new=NULL}
  }


  if(method == "gpr"){
    ## Fitting P(-C<=-c_0|X) (since P(C>=c_0|X)=P(-C<=-c_0|X))
    fit$C <- -fit$C
    gpr_mdl <- GauPro(X = as.matrix(fit[,names(fit) %in% xnames]),
                    Z = fit$C, D = p,
                    type = "Gauss")
    
    ## Computing the censoring scores for the fitting data
    mean_fit <- gpr_mdl$predict(fit[,names(fit) %in% xnames])
    sd_fit <- gpr_mdl$predict(fit[,names(fit) %in% xnames],
                              se.fit = TRUE)$se

    pr_fit <- pnorm((-c - mean_fit) / sd_fit)

    ## Computing the censoring scores for the calibration data
    mean_calib <- gpr_mdl$predict(calib[,names(calib) %in% xnames])
    sd_calib <- gpr_mdl$predict(calib[,names(calib) %in% xnames],
                              se.fit = TRUE)$se

    pr_calib <- pnorm((-c - mean_calib) / sd_calib)

    ## Computing the censoring scores for the test data
    if(!is.null(test)){
      newdata <- data.frame(test)
      colnames(newdata) <- xnames
      mean_new <- gpr_mdl$predict(newdata[,names(newdata) %in% xnames])
      sd_new <- gpr_mdl$predict(newdata[,names(newdata) %in% xnames],
                              se.fit = TRUE)$se

      pr_new <- pnorm((-c - mean_new) / sd_new)

    }else{pr_new=NULL}

  }
  return(list(pr_fit = pr_fit, pr_calib = pr_calib, pr_new = pr_new))
}


#' Computing the calibration term with covaraite shift
#'
#' construct the one-sided confidence interval for a unit's survival time T
#'
#' @param x a vector of the covariate of the test data.
#' @param r the censoring time of the test data.
#' @param alpha a number betweeo 0 and 1, specifying the miscaverage rate.
#' @param data a data frame used for calibration, containing four columns: (X,R,event,censored_T). 
#' @param mdl The fitted model to estimate the conditional quantile (default is NULL).
#' @param quant_lo the fitted conditional quantile for the calibration data (default is NULL).
#' @param new_quant_lo the fitted conditional quantile for the test data (default is NULL).
#'
#' @return low_ci a value of the lower bound for the survival time of the test point.
#' @return includeR 0 or 1, indicating if [r,inf) is included in the confidence interval.
#'
#' @family confint
#'
#' @export


get_calibration <- function(score,weight_calib,weight_new,alpha){
  ## Check input format
  if(length(score)!=length(weight_calib)) stop("The length of score is not compatible with the length of weight!")

  if(!is.numeric(alpha)) stop("alpha should be a real number between 0 and 1!")
  if(alpha>1 | alpha<0) stop("alpha should be a real number between 0 and 1!")

  ## Computing the calibration term
  weight <- c(weight_calib,weight_new)
  weight <- weight/sum(weight)
  score_vec <- c(score,Inf)
  sort_score <- sort(score_vec)
  order_score <- order(score_vec)
  sort_w <- weight[order_score]
  idxw <- min(which(cumsum(sort_w)>=1-alpha))
  calib_term <- sort_score[idxw]

  return(calib_term)
}

 
