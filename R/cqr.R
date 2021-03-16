#' Conformalized Quantile Regression
#'
#' @export

cqr <- function(x,Xtrain,Ytrain,
                alpha=0.05,
                I_fit=NULL,seed=24601){

  X <- Xtrain
  Y <- Ytrain
  set.seed(seed)
  if(is.null(dim(X)[1])){
    n <- length(X)
    p <- 1
  }else{
    n <- dim(X)[1]
    p <- dim(X)[2]
  }
  xnames <- paste0('X', 1:p)
  data <- as.data.frame(cbind(Y,X))
  colnames(data) <- c("Y",xnames)

  ## Divide the dataset for model fitting and calibration 
  n_train = n/2
  n_calib = n-n_train
  if(is.null(I_fit)){
    I_fit <- sample(1:n,n_train,replace = FALSE)
  }
  data_fit <- data[I_fit,]
  data_calib <- data[-I_fit,]
  newdata <- data.frame(x)
  colnames(newdata) <- xnames

  ## Run quantile regression on data_fit
  fmla <- as.formula(paste("Y~ ", paste(xnames, collapse= "+")))
  mdl <- rq(fmla,data = data_fit,tau = alpha)

  ## Obtain the scores
  res <- predict(mdl,
                 newdata = data_calib,
                 type="quantile")
  quant <-  res  
  score <- quant-data_calib$Y

  ## Get the fitted quantile for the test point 
  res <- predict(mdl,
                 newdata = newdata,
                 type="quantile")
  new_quant <-  res  
  
  ## Get the calibration term
  corr_term <- quantile(c(score,Inf),1-alpha)
  lower_bnd <- new_quant-corr_term
  return(lower_bnd)


}
