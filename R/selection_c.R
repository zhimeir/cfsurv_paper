#' Automatic selection of censoring time
#'
#' The funciton that automatically choose a value of c_0.
#'
#' @export

selection_c <- function(X,C,event,time,alpha,
                        c_ref,weight_ref,
                        model="cox",
                        type="quantile",
                        dist="weibull"){
  
  ## Get the dimension of the input
  if(is.null(dim(X))){
    n <- length(X)
    p <- 1
  }else{
    n <- dim(X)[1]
    p <- dim(X)[2]
  }
  xnames <- paste0("X",1:p)
  data <- cbind(X,C,event,time)
  data <- data.frame(data)
  colnames(data) <- c(xnames,"C","event","censored_T")

  ## Evaluate the average bound for each candidate c
  bnd_ref <- c()
  for(i in 1:length(c_ref)){
  bnd <- evaluate_length(c_ref[i],alpha=alpha,n=n,p=p,model,
                        data=data,weight=weight_ref[,i],xnames=xnames,
                        type=type,dist=dist)
  bnd_ref <- c(bnd_ref,bnd)
  }
  c_opt <- c_ref[which.max(bnd_ref)]
  return(list(c_opt=c_opt,c_ref=c_ref,bnd_ref=bnd_ref))

}

evaluate_length <- function(c,alpha,n,p,
                            model,
                            data,
                            weight,
                            xnames,
                            type = "quantile",
                            dist = "weibull",
                            seed = 2020){
  ## Determine the fitting set, the calibration set and the test set
  set.seed(seed)
  I_fit <- sample(1:n,floor(n/2),replace=FALSE)
  I_calib <- sample((1:n)[-I_fit],floor(n/4),replace=FALSE)
  I_test <- (1:n)[-c(I_fit,I_calib)]
  
  data_fit <- data[I_fit,]
  data_calib <- data[I_calib,]
  data_test <- data[I_test,]
  
  if(is.null(weight)){
    res <- censoring_prob(fit=data_fit,calib=data_calib,test=data_test,
                          xnames=xnames,c)
    pr_calib <- res$pr_calib
    pr_new <- res$pr_new
    weight_calib <- 1/pr_calib
    weight_new <- 1/pr_new
  }else{
    weight_calib <- weight[I_calib]
    weight_new <- weight[I_test]
  }
  x <- data_test[,colnames(data_test)%in%xnames]
  
  if(model == "cox"){
    bnd <- cox_based(x,c,alpha,
                     data_fit,
                     data_calib,
                     type = "quantile",
                     dist,
                     weight_calib,
                     weight_new)
   }
  
  if(model == "randomforest"){
    bnd <- rf_based(x,c,alpha,
                    data_fit,
                    data_calib,
                    weight_calib,
                    weight_new)
  }
  
  if(model == "pow"){
    bnd <- pow_based(x,c,alpha,
                    data_fit,
                    data_calib,
                    weight_calib,
                    weight_new)
  }

  if(model == "portnoy"){
    bnd <- portnoy_based(x,c,alpha,
                        data_fit,
                        data_calib,
                        weight_calib,
                        weight_new)
  }

  if(model == "PengHuang"){
    bnd <- ph_based(x,c,alpha,
                   data_fit,
                   data_calib,
                   weight_calib,
                   weight_new)
  }

  return(mean(bnd))
}
