#!/usr/bin/env Rscript
########################################
## Start of problem independent section
########################################
args <- commandArgs(trailingOnly = TRUE)
seed<- as.integer(args[1])
if(is.na(seed)){seed <- 1}
cat(sprintf(" - Running the script with seed %d.\n", seed))

########################################
## load libraries
########################################
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(quantreg))
suppressPackageStartupMessages(library(conTree))
suppressPackageStartupMessages(library(GauPro))
suppressPackageStartupMessages(library(gbm))

########################################
### source code
########################################
source("../utils/source.R")

########################################
## Output directory
########################################
out_dir <- "../results/ld_homosc"
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}

########################################
## Parameter
########################################
n <- 3000
n_train <- n / 2
n_calib <- n / 2
n_test <- 3000
beta <- 20 / sqrt(n)
c_ref <- 1 : 6 / 2
xmin <- 0; xmax <- 4
exp_rate <- 0.4
pr_all_list <- matrix(0, n+n_test, length(c_ref))
alpha_list <- 0.1

########################################
## Data generating models
########################################
gen_t <- function(x) exp(2 + beta * sqrt(abs(x)) +  1.5 * rnorm(length(x))) 
gen_c <- function(x) rexp(rate = exp_rate, n = length(x)) 

########################################
## Generate training data
########################################
set.seed(24601)
X <- runif(n_train, xmin, xmax)
T <- gen_t(X)
C <- gen_c(X)
event <- (T < C)
censored_T <- pmin(T, C)
data_fit <- data.frame(X1 = X, C = C, censored_T = censored_T, event = event)

########################################
## Generate the calibration data and the test data
########################################
set.seed(seed)
X <- runif(n_calib + n_test, xmin, xmax)     
T <- gen_t(X) 
C <- gen_c(X)
event <- (T < C)
censored_T <- pmin(T, C)
data <- data.frame(X1 = X, C = C, event = event, censored_T = censored_T)
data_calib <- data[1 : n_calib, ]
data_test <- data[(n_calib + 1) : (n_calib + n_test), ]
data <- rbind(data_fit, data_calib)

########################################
## Estimate pr(C>=c)
########################################
xnames <- paste0("X",1) 
for(i in 1:length(c_ref)){
  res <- censoring_prob(fit = data_fit, calib = data_calib,
                        test = data_test[,colnames(data_test)%in%xnames],
                        xnames = xnames, c = c_ref[i], 
                        method = "distBoost", n.tree = 500)
  pr_all_list[,i] <- c(res$pr_fit, res$pr_calib, res$pr_new)
}
pr_list <- pr_all_list[1 : n, ]
pr_new_list <- pr_all_list[(n + 1) : (n + n_test), ]

########################################
## CDR + conTree
########################################
cat(" - Computing the result of CDR-LPB - distBoost...")
res <- lapply(alpha_list,
              cfsurv,
              x=data_test$X1,
              c_list = c_ref,
              Xtrain = data$X1,
              C = data$C,
              event = data$event,
              time=data$censored_T,
              I_fit = 1:n_train,
              pr_list = pr_list,
              pr_new_list = pr_new_list,
              model="distBoost",
              seed = seed+7,
              n.tree= 500
)
res <- do.call(rbind,lapply(res,as.data.frame))
output <- data.frame(distboost.bnd = res[,1])
cat("done.\n")

########################################
## CQR + conTree
########################################
cat(" - Computing the result of CQR-LPB - conTree...")
res <- lapply(alpha_list,
              cfsurv,
              x=data_test$X1,
              c_list = c_ref,
              Xtrain = data$X1,
              C = data$C,
              event = data$event,
              time=data$censored_T,
              I_fit = 1:n_train,
              pr_list = pr_list,
              pr_new_list = pr_new_list,
              model="quantBoost",
              seed = seed+7,
              n.tree = 500
)
res <- do.call(rbind,lapply(res,as.data.frame))
output$quantboost.bnd <- res[,1]
cat("done.\n")

########################################
## vanilla CQR
########################################
cat(" - Computing the result of vanilla-CQR...")
res <- lapply(alpha_list,
              cqr,
              x = data_test$X1,
              Xtrain = data$X1,
              Ytrain = data$censored_T,
              I_fit = 1 : n_train,
              seed = seed + 7)
res <- do.call(rbind,lapply(res,as.data.frame))
output$cqr.bnd <- res[,1]
cat("done.\n")


########################################
## CQR + censored random forest
########################################
cat(" - Computing the result of CQR-LPB - randomForest...")
res <- lapply(alpha_list,
              cfsurv,
              x=data_test$X1,
              c_list= c_ref,
              Xtrain = data$X1,
              C = data$C,
              event = data$event,
              time=data$censored_T,
              I_fit = 1:n_train,
              type="quantile",
              pr_list = pr_list,
              pr_new_list = pr_new_list,
              model="randomforest",
              seed = seed+7,
              n.tree = 500
              )
res <- do.call(rbind,lapply(res,as.data.frame))
output$cqr.rf.bnd <- res[,1]
cat("done.\n")


########################################
### Utility functions
########################################
## A function to get result
extract_res_univariate <- function(x){
  res <- mdl_coef[1] + mdl_coef[-1]*x
  return(res)
}

## A utility function to extract quantiles from a coxph object
extract_quant <- function(mdl, x, alpha){
  res <- summary(survfit(mdl, newdata = x))
  time_point <- res$time
  survcdf <- 1 - res$surv
  quant <- time_point[min(which(survcdf >= alpha))]
  return(quant)
}
########################################
## Quantile regression
########################################

## Cox
cat(" - Computing the result of Cox...")
mdl <- coxph(Surv(censored_T,event)~X1,data = data)
res <- c()
for(alpha in alpha_list){
  res <- c(res, apply(data_test, 1, extract_quant, mdl = mdl,  alpha = alpha))
}
output$cox.bnd <- res
cat("Done.\n")

## AFT + Weibull
cat(" - Computing the result of Cox...")
mdl <- survreg(Surv(censored_T,event)~X1,data = data, dist = "weibull")
res <- predict(mdl,
               newdata = data_test,
               type = "quantile",
               p = alpha_list)
output$aft.bnd <- as.vector(res)
cat("Done.\n")

## Powell
cat(" - Computing the result of powell...")
mdl <- crq(Curv(censored_T,C,ctype = "right")~X1,data=data,taus = alpha_list,method = "Pow")
res <- predict(mdl,
        newdata = data_test,
        type = "quantile")
output$pow.bnd <- as.vector(res)
cat("Done.\n")

## Portnoy 
cat(" - Computing the result of portnoy...")
mdl <- crq(Surv(censored_T,event)~X1,data = data, method = "Portnoy")
mdl_coef <- coef(mdl,taus = alpha_list)
res <- sapply(data_test$X1,extract_res_univariate)
output$por.bnd <- as.vector(res)
cat("Done.\n")

## Peng and Huang
cat(" - Computing the result of Peng and Huang...")
mdl <- crq(Surv(censored_T,event)~X1,data = data, method = "PengHuang")
mdl_coef <- coef(mdl,taus = alpha_list)
res <- sapply(data_test$X1,extract_res_univariate)%>%t()
output$penghuang.bnd <- as.vector(res)
cat("Done.\n")

## random forest
cat(" - Computing the result of random forest...")
ntree <- 1000
nodesize <- 80
res <- c()
for(alpha in alpha_list){
  mdl <- crf.km(as.formula("censorerd_T~X1"),
              ntree = ntree, 
              nodesize = nodesize,
              data_train = data[,colnames(data) %in% c("X1","censored_T","event")], 
              data_test = data.frame(X1 = data_test$X1), 
              yname = 'censored_T', 
              iname = 'event',
              tau = alpha,
              method = "grf")
  res <- c(res,mdl$predicted)
}
output$rf.bnd  <- res
cat("Done.\n")

########################################
## Append the ground truth
########################################
output$T <- rep(T[(n_calib + 1) : (n_calib + n_test)],length(alpha_list))
output$R <- rep(data_test$C,length(alpha_list))
output$X <- rep(data_test$X1,length(alpha_list))

########################################
## Save the result
########################################
save_dir <- paste0(out_dir,"res",seed,".csv")
write_csv(output,save_dir)
