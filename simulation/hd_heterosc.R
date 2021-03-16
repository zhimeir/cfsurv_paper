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
## Output direcroty
########################################
out_dir <- "../results/hd_heterosc"
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
p <- 100
beta <- 30 / sqrt(n)
xnames <- paste0("X",1:p) 
c_ref <-1 : 6 / 2
exp_rate <- 0.4
alpha_list <- 0.1

########################################
## Data generating models
########################################
mu_x <- function(x) beta * x[,1]^2 - beta * x[,3] * x[,5] + 1
sigma_x <- function(x) (abs(x[,10]) + 1) 
gen_t <- function(x) 2 * exp(mu_x(x) + sigma_x(x) * rnorm(dim(x)[1]))
gen_c <- function(x) rexp(rate = exp_rate, n = dim(x)[1])

## Generate training data
set.seed(24601)
X <- matrix(runif(n_train * p, min = -1, max = 1), n_train)
T <- gen_t(X)
C <- gen_c(X) 
event <- (T<C)
censored_T <- pmin(T,C)
data_fit <- data.frame(X, C = C, censored_T = censored_T, event = event)
colnames(data_fit) <- c(xnames, "C", "censored_T", "event")

########################################
## Generate the calibration data and the test data
########################################
set.seed(seed)
X <- matrix(runif((n_calib + n_test) * p, min = -1, max = 1), n_calib + n_test)
T <- gen_t(X)
C <- gen_c(X)
event <- (T<C)
censored_T <- pmin(T,C)
data <- data.frame(X, C = C, censored_T = censored_T,  event = event)
colnames(data) <- c(xnames, "C", "censored_T", "event")
data_calib <- data[1:n_calib,]
data_test <- data[(n_calib+1) : (n_calib+n_test),]
data <- rbind(data_fit,data_calib)

########################################
## Estimate pr(C>=c)
########################################
pr_all_list <- matrix(0, n + n_test, length(c_ref))
for(i in 1 : length(c_ref)){
  res <- censoring_prob(fit = data_fit, calib = data_calib,
                        test = data_test[,colnames(data_test)%in%xnames],
                        xnames = xnames, c = c_ref[i], method = "distBoost", n.tree = 500)
  pr_all_list[,i] <- c(res$pr_fit, res$pr_calib, res$pr_new)
}
pr_list <- pr_all_list[1 : n, ]
pr_new_list <- pr_all_list[(n + 1) : (n + n_test), ]


########################################
### CDR + conTree 
########################################
cat("Computing the results of CDR + conTree...")
res <- lapply(alpha_list,
              cfsurv,
              x=data_test[,colnames(data_test)%in%xnames],
              c_list = c_ref,
              Xtrain = data[,colnames(data)%in%xnames],
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
cat(" - Computing the results of CQR + conTree...")
res <- lapply(alpha_list,
              cfsurv,
              x = data_test[,colnames(data_test)%in%xnames],
              c_list = c_ref,
              Xtrain = data[,colnames(data)%in%xnames],
              C = data$C,
              event = data$event,
              time=data$censored_T,
              I_fit = 1:n_train,
              pr_list = pr_list,
              pr_new_list = pr_new_list,
              model="quantBoost",
              seed = seed+7,
              n.tree= 500
)
res <- do.call(rbind,lapply(res,as.data.frame))
output$quantboost.bnd <- res[,1]
cat("done.\n")

########################################
## vanilla CQR
########################################
cat("Computing the result of vanilla CQR...")
res <- lapply(alpha_list,
              cqr,
              x=data_test[,colnames(data_test)%in%xnames],
              Xtrain = data[,colnames(data)%in%xnames],
              Ytrain = data$censored_T,
              I_fit = 1:n_train,
              seed = seed +7)
res <- do.call(rbind,lapply(res,as.data.frame))
output$cqr.bnd <- res[,1]
cat("done.\n")

########################################
## CQR + censored random forest
########################################
cat("Computing the result of censoredCQR + randomforest...")
res <- lapply(alpha_list,
              cfsurv,
              x=data_test[,colnames(data_test)%in%xnames],
              c_list= c_ref,
              Xtrain = data[,colnames(data)%in%xnames],
              C = data$C,
              event = data$event,
              time=data$censored_T,
              I_fit = 1:n_train,
              type="quantile",
              pr_list = pr_list,
              pr_new_list = pr_new_list,
              model="randomforest",
              seed = seed+7
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
fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
            paste(xnames, collapse= "+")))
mdl <- coxph(fmla, data = data)
res <- c()
for(alpha in alpha_list){
  res <- c(res, apply(data_test, 1, extract_quant, mdl = mdl,  alpha = alpha))
}
output$cox.bnd <- res
cat("done.\n")

## AFT + Weibull
cat(" - Computing the result of Cox...")
fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
            paste(xnames, collapse= "+")))
mdl <- survreg(fmla, data = data, dist = "weibull")
res <- predict(mdl,
               newdata = data_test,
               type = "quantile",
               p = alpha_list)
output$aft.bnd <- as.vector(res)
cat("done.\n")

## Powell
cat(" - Computing the result of powell...")
fmla <- as.formula(paste("Curv(censored_T, C, ctype = \"right\") ~ ",
            paste(xnames, collapse= "+")))
mdl <- crq(fmla,data=data,taus = alpha_list,method = "Pow")
res <- predict(mdl,
        newdata = data_test,
        type = "quantile")
output$pow.bnd <- as.vector(res)
cat("done.\n")

## Portnoy 
cat(" - Computing the result of portnoy...")
fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
            paste(xnames, collapse= "+")))
mdl <- crq(fmla,data = data, method = "Portnoy")
mdl_coef <- coef(mdl,taus = alpha_list)
res <- as.matrix(data_test[, colnames(data_test) %in% xnames]) %*% 
  mdl_coef[-1] + mdl_coef[1] 
output$por.bnd <- as.vector(res)
cat("done.\n")

## Peng and Huang
cat(" - Computing the result of Peng and Huang...")
fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
            paste(xnames, collapse= "+")))
mdl <- crq(fmla,data = data, method = "PengHuang")
mdl_coef <- coef(mdl,taus = alpha_list)
res <- as.matrix(data_test[, colnames(data_test) %in% xnames]) %*% 
  mdl_coef[-1] + mdl_coef[1] 
output$penghuang.bnd <- as.vector(res)
cat("done.\n")

## random forest
cat(" - Computing the result of random forest...")
ntree <- 1000
nodesize <- 80
fmla <- as.formula(paste("censored_T~ ",
            paste(xnames, collapse= "+")))
res <- c()
for(alpha in alpha_list){
  mdl <- crf.km(fmla,
              ntree = ntree, 
              nodesize = nodesize,
              data_train = data[,colnames(data) %in% c(xnames,"censored_T","event")], 
              data_test = data_test[,colnames(data_test) %in% xnames], 
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
output$T <- T[(n_calib + 1) : (n_calib + n_test)]
output$C <- data_test$C
output <- cbind(output, data_test[,colnames(data_test) %in% paste0("X",c(1,2,3,5,10))])

########################################
## Save the result
########################################
save_dir <- paste0(out_dir,"res",seed,".csv")
write_csv(output,save_dir)
