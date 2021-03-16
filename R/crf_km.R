#' crf.km
#'
#' The function is based on the scripts from https://github.com/AlexanderYogurt/censored_ExtremelyRandomForest to perform censored quantile random forest. 
#'
#' @export


crf.km <- function(fmla, ntree, 
                   nodesize,
                   data_train,
                   data_test,
                   yname, 
                   iname,
                   tau,
                   xname = NULL, 
                   calibrate_taus = c(0.1, 0.5, 0.9),
                   honesty = TRUE, 
                   method = "grf", 
                   splitrule = "extratrees", 
                   nnb = FALSE, reg.split = FALSE) {
  # build Forest model
  if (method == "randomForest") {
    rf <- randomForest(fmla, data = data_train, ntree = ntree, nodesize=nodesize)
    # get proximity matrix
    proxMtx <- rf.getWeights(rf, data_train, data_test, yname, iname)
  } else if (method == "ranger") {
    if (is.null(xname)) {
      rf <- ranger(fmla, data = data_train, num.trees = ntree, min.node.size=nodesize, splitrule=splitrule)
    } else {
      rf <- ranger(fmla, data = data_train[,c(xname,yname)], num.trees = ntree, min.node.size=nodesize, splitrule=splitrule)
    }
    # get proximity matrix
    proxMtx <- ranger.getWeights(rf, data_train, data_test, yname, iname, xname)
  } else if (method == "grf") {
    if (is.null(xname)) {
      rf <- quantile_forest(X=data_train[ ,!(names(data_train) %in% c(yname, iname)), drop=F],
                            Y=data_train[,yname],
                            quantiles = calibrate_taus, 
                            num.trees = ntree,
                            min.node.size = nodesize,
                            honesty = honesty,
                            regression.splitting = reg.split)
    } else {
      rf <- quantile_forest(X=data_train[ ,xname], Y=data_train[,yname], quantiles = calibrate_taus, 
                            num.trees = ntree, min.node.size = nodesize, honesty = honesty,
                            regression.splitting = reg.split)
    }
    proxMtx <- as.matrix(grf.getWeights(rf, data_test, yname, iname, xname))
  }
  
  # censor forest
  n <- nrow(data_test)
  ntrain <- nrow(data_train)
  Yc <- rep(NA, n) # to store new predictions
  Ytrain <- data_train[[yname]]
  censorInd <- data_train[[iname]]
  right_censoring <- TRUE
  
  # find minimum
  for (r in 1:NROW(data_test)) {
    # C survival estimate
    if (nnb) {
      boot.idx <- proxMtx[r, ] > 1/ntrain
      C.boot <- Ytrain[boot.idx]
      i.boot <- 1 - censorInd[boot.idx]
      C.km <- survfit(Surv(C.boot, i.boot) ~ 1, type = 'kaplan-meier')
      C.surv <- stepfun(C.km$time, c(1, C.km$surv))
    } else {
      boot.idx <- proxMtx[r, ] > 0
    }

    max.uncensored <- max(Ytrain[boot.idx & censorInd==1])
    candidateY <- sort(Ytrain[boot.idx & Ytrain<=max.uncensored], decreasing = FALSE)

    Yc[r] <- candidateY[1]
    min_loss <- 1000000
    denom <- sapply(Ytrain, function(x) {1*(Ytrain >= x)%*%proxMtx[r, ]})
    denom[denom==0] <- 1
    base <- 1 - (proxMtx[r, ]/denom)
    base <- base^(1-censorInd)
    nn <- 1
    for (lambda in candidateY) {
      if (right_censoring) {
        if (nnb) {
          kappa <- C.surv(lambda)
        } else {
          kappa <- C.surv(lambda, Ytrain, base)
          #print(kappa)
        }
        loss1 <- proxMtx[r, ]%*%(1*(Ytrain > lambda))
        #print(kappa)
        loss <- abs((1-tau)*kappa - loss1)
      } else {
        #
      }
      if (loss < min_loss) {
        Yc[r] <- lambda
        min_loss <- loss
      } else if (loss == min_loss) {
        Yc[r] <- (Yc[r]*nn + lambda) / (nn+1)
        nn <- nn + 1
      }
    }
  }
  return(
    list(
      'predicted'=Yc,
      'proxMtx'=proxMtx
    )
  )
}

#' C.surv
#'
#' An internal function for censored quantile random forest
#'
#' @export
C.surv <- function(q, Y, base) {
  result <- prod(base[Y <= q])
  return(result)
}


#' grf.getWeights
#'
#' An internal function for censored quantile random forest
#'
#' @export
grf.getWeights = function(grf, testdata, y_name, c_name, x_name=NULL) {
  
  if (is.null(x_name)){
    weights <- get_sample_weights(grf, newdata=testdata[ ,!(names(testdata) %in% c(y_name, c_name)), drop=F]) 
  } else {
    weights <- get_sample_weights(grf, newdata=testdata[ ,x_name])
  }
  
  return (weights)
}


#' rf.getNodes
#'
#' An internal function for censored quantile random forest
#'
#' @export
rf.getNodes = function(rf, data, y_name, c_name) {
  pred = predict(rf, data[ ,!(names(data) %in% c(y_name, c_name)), drop=F], nodes = T)
  nodes = attributes(pred)$nodes

  return(nodes)
}


#' rf.getWeights
#'
#' An internal function for censored quantile random forest
#'
#' @export
rf.getWeights = function(rf, traindata, testdata, y_name, c_name) {
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  cat("number of cores: ", cores[1]-1)
  registerDoParallel(cl)

  # retrieve training nodes
  nodes = rf.getNodes(rf, traindata, y_name, c_name) # [train.samples, trees]

  # retrieve nodes for test data
  test.nodes = rf.getNodes(rf, testdata, y_name, c_name) # [test.samples, trees]
  nodesize = test.nodes # [test.samples, trees]

  # loop over trees
  print("loop begins...\n")
  ntrees = rf$ntree
  weights <- foreach(k=1:ntrees, .combine='+') %dopar% {
      in.leaf = 1*outer(test.nodes[,k],nodes[,k],'==') # [test.samples, train.samples] indicates whether a test and a train data are in the same leaf node
      mapping=plyr::count(nodes[,k])
      nodesize[,k] = plyr::mapvalues(test.nodes[,k], from=mapping$x, to=mapping$freq, warn_missing=F)
      weight = in.leaf / nodesize[,k]
      weight
  }
  stopCluster(cl)

  return(weights/ntrees)
}


#' ranger.getNodes
#'
#' An internal function for censored quantile random forest
#'
#' @export
ranger.getNodes = function(rg, data, y_name, c_name, x_name=NULL) {
  if (is.null(x_name)) {
    pred = predict(rg, data[ ,!(names(data) %in% c(y_name, c_name)), drop=F], type = "terminalNodes")
  } else {
    pred = predict(rg, data[ ,x_name], type = "terminalNodes")
  }
  nodes = pred$predictions

  return(nodes)
}


#' ranger.getWeights
#'
#' An internal function for censored quantile random forest
#'
#' @export
ranger.getWeights = function(rg, traindata, testdata, y_name, c_name, x_name=NULL) {
  num_cores = detectCores() - 1 #not to overload your computer
  cl <- makeCluster(num_cores, type="FORK")
  cat("number of cores: ", num_cores)
  registerDoParallel(cl)

  # retrieve training nodes
  #print("retrieving training node information...\n")
  nodes = ranger.getNodes(rg, traindata, y_name, c_name, x_name) # [train.samples, trees]

  # retrieve nodes for test data
  #print("retrieving test node information...\n")
  test.nodes = ranger.getNodes(rg, testdata, y_name, c_name, x_name) # [test.samples, trees]
  nodesize = test.nodes # [test.samples, trees]

  # loop over trees
  print("loop begins...")
  ntrees = rg$num.trees
  weights <- foreach(k=1:ntrees, .combine='+') %dopar% {
    in.leaf = 1*outer(test.nodes[,k],nodes[,k],'==') # [test.samples, train.samples] indicates whether a test and a train data are in the same leaf node
    mapping = plyr::count(nodes[,k])
    nodesize[,k] = plyr::mapvalues(test.nodes[,k], from=mapping$x, to=mapping$freq, warn_missing=F)
    weight = in.leaf / nodesize[,k]
    weight
  }
  
  stopCluster(cl)
  stopImplicitCluster()

  return(weights/ntrees)
}
