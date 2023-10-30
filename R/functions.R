
##### transreg #####

#' @export
#' 
#' @title
#' Penalised regression with multiple sets of prior effects
#' 
#' @description
#' Implements penalised regression with multiple sets of prior effects
#' 
#' @param y
#' target: vector of length \eqn{n} (see \code{family})
#' @param X
#' features: matrix with \eqn{n} rows (samples)
#' and \eqn{p} columns (features)
#' @param prior
#' prior coefficients: matrix with \eqn{p} rows (features)
#' and \eqn{k} columns (sources of co-data)
#' @param family
#' character "gaussian" (\eqn{y}: real numbers),
#' "binomial" (\eqn{y}: 0s and 1s),
#' or "poisson" (\eqn{y}: non-negative integers);
#' @param alpha
#' elastic net mixing parameter (0=ridge, 1=lasso):
#' number between 0 and 1
#' @param foldid
#' fold identifiers: vector of length \eqn{n}
#' with entries from 1 to `nfolds`
#' @param nfolds
#' number of folds: positive integer
#' @param scale
#' character
#' "exp" for exponential calibration or
#' "iso" for isotonic calibration
#' @param stack
#' character "sta" (standard stacking) or "sim" (simultaneous stacking)
#' @param sign
#' sign discovery procedure: logical
#' (experimental argument)
#' @param switch
#' choose between positive and negative weights for each source: logical
#' @param select
#' select from sources: logical
#' @param track
#' show intermediate output (messages and plots): logical
#' @param parallel
#' logical (see cv.glmnet)
#'
#' @details
#' * \eqn{n}: sample size
#' * \eqn{p}: number of features
#' * \eqn{k}: number of sources
#' 
#' @return
#' Returns an object of class `transreg`.
#' Rather than accessing its slots (see list below),
#' it is recommended to use methods like
#' [coef.transreg()] and [predict.transreg()].
#' 
#' * slot `base`:
#' Object of class `glmnet`.
#' Regression of outcome on features (without prior effects),
#' with \eqn{1 + p} estimated coefficients
#' (intercept + features).
#' 
#' * slot `meta.sta`:
#' `NULL` or object of class `glmnet`.
#' Regression of outcome on cross-validated linear predictors
#' from prior effects and estimated effects,
#' with \eqn{1 + k + 2} estimated coefficients
#' (intercept + sources of co-data + lambda_min and lambda_1se).
#' 
#' * slot `meta.sim`:
#' `NULL` or object of class `glmnet`.
#' Regression of outcome on meta-features
#' (cross-validated linear predictors from prior effects)
#' and original features,
#' with \eqn{1 + k + p} estimated coefficients
#' (intercept + sources of co-data + features).
#' 
#' * slot `prior.calib`:
#' Calibrated prior effects.
#' Matrix with \eqn{p} rows and \eqn{k} columns.
#' 
#' * slot `data`:
#' Original data.
#' List with slots `y`, `X` and `prior` (see arguments).
#' 
#' * slot `info`:
#' Information on call.
#' Data frame with entries
#' \eqn{n}, \eqn{p}, \eqn{k}, `family`, `alpha`, `scale` and `stack`
#' (see details and arguments).
#' 
#' @inherit transreg-package references
#' 
#' @seealso
#' Methods for objects of class [transreg]
#' include \code{\link[=coef.transreg]{coef}} 
#' and \code{\link[=predict.transreg]{predict}}.
#' 
#' @examples
#' #--- simulation ---
#' n <- 100; p <- 500
#' X <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
#' beta <- rnorm(p)*rbinom(n=p,size=1,prob=0.2)
#' prior1 <- beta + rnorm(p)
#' prior2 <- beta + rnorm(p)
#' y_lin <- X %*% beta
#' y_log <- 1*(y_lin > 0)
#' 
#' #--- single vs multiple priors ---
#' one <- transreg(y=y_lin,X=X,prior=prior1)
#' two <- transreg(y=y_lin,X=X,prior=cbind(prior1,prior2))
#' weights(one)
#' weights(two)
#' 
#' \donttest{
#' #--- linear vs logistic regression ---
#' lin <- transreg(y=y_lin,X=X,prior=prior1,family="gaussian")
#' log <- transreg(y=y_log,X=X,prior=prior1,family="binomial")
#' hist(predict(lin,newx=X)) # predicted values
#' hist(predict(log,newx=X)) # predicted probabilities
#' 
#' #--- ridge vs lasso penalisation ---
#' ridge <- transreg(y=y_lin,X=X,prior=prior1,alpha=0)
#' lasso <- transreg(y=y_lin,X=X,prior=prior1,alpha=1)
#' # initial coefficients (without prior)
#' plot(x=coef(ridge$base)[-1]) # dense
#' plot(x=coef(lasso$base)[-1]) # sparse
#' # final coefficients (with prior)
#' plot(x=coef(ridge)$beta) # dense
#' plot(x=coef(lasso)$beta) # not sparse
#' 
#' #--- exponential vs isotonic calibration ---
#' exp <- transreg(y=y_lin,X=X,prior=prior1,scale="exp")
#' iso <- transreg(y=y_lin,X=X,prior=prior1,scale="iso")
#' plot(x=prior1,y=exp$prior.calib)
#' plot(x=prior1,y=iso$prior.calib)
#' 
#' #--- standard vs simultaneous stacking ---
#' prior <- c(prior1[1:250],rep(0,250))
#' sta <- transreg(y=y_lin,X=X,prior=prior,stack="sta")
#' sim <- transreg(y=y_lin,X=X,prior=prior,stack="sim")
#' plot(x=coef(sta$base)[-1],y=coef(sta)$beta)
#' plot(x=coef(sim$base)[-1],y=coef(sim)$beta)}
#' 
transreg <- function(y,X,prior,family="gaussian",alpha=1,foldid=NULL,nfolds=10,scale="iso",stack="sim",sign=FALSE,switch=FALSE,select=TRUE,track=FALSE,parallel=FALSE){
  
  if(is.vector(prior)){
    prior <- matrix(prior,ncol=1)
  }
  prior[is.na(prior)] <- 0 # allows for missing prior coefficients
  if(all(y %in% c(0,1)) != (family=="binomial")){stop("Check 'family' of 'y'.")}
  if(length(y)!=nrow(X)){stop("Entries in 'y' must match rows in 'X'.")}
  if(ncol(X)!=nrow(prior)){stop("Columns in 'X' must match rows in 'prior'.")}
  if(all(prior>=0) & !sign){warning("Consider sign discovery procedure.")}

  foldid <- .folds(y=y,nfolds=nfolds,foldid=foldid)$foldid

  n <- nrow(X); p <- ncol(X)
  k <- ncol(prior)
  
  base <- glmnet::cv.glmnet(y=y,x=X,family=family,alpha=alpha,nlambda=100,keep=TRUE,foldid=foldid,parallel=parallel)
  
  # external re-scaling
  prior.ext <- prior
  if(sign){
    prior.ext <- .signdisc(y=y,X=X,prior=prior.ext,family=family)
  }
  base$z <- prior.ext
  if(scale=="exp"){
    prior.ext <- .exp.multiple(y=y,X=X,prior=prior.ext,family=family,switch=switch,select=select,track=track)
  } else if(scale=="iso"){
    prior.ext <- .iso.multiple(y=y,X=X,prior=prior.ext,family=family,switch=switch,select=select,track=track)
  } else {
    stop("Invalid scale.",call.=FALSE)
  }
  
  y_hat <- matrix(NA,nrow=n,ncol=k+2)

  for(i in seq_len(nfolds)){
    y0 <- y[foldid!=i]
    X0 <- X[foldid!=i,]
    X1 <- X[foldid==i,]
    
    # internal re-scaling
    prior.int <- prior
    if(sign){
      prior.int <- .signdisc(y=y0,X=X0,prior=prior.int,family=family)
    }
    if(scale=="exp"){
      prior.int <- .exp.multiple(y=y0,X=X0,prior=prior.int,family=family,switch=switch,select=select,track=track)
    } else if(scale=="iso"){
      prior.int <- .iso.multiple(y=y0,X=X0,prior=prior.int,family=family,switch=switch,select=select,track=track)
    } else {
      stop("Invalid scale.",call.=FALSE)
    }
    
    # predictions
    for(j in seq_len(k)){
      y_hat[foldid==i,j] <- X1 %*% prior.int$beta[,j]
    }
    y_hat[foldid==i,k+1] <- base$fit.preval[base$foldid==i,base$lambda==base$lambda.min]
    y_hat[foldid==i,k+2] <- base$fit.preval[base$foldid==i,base$lambda==base$lambda.1se]
  }
  
  base$prior <- prior.ext
  
  #--- standard stacking ---
  if(any(stack=="sta")){
     meta.sta <- glmnet::cv.glmnet(y=y,x=y_hat,foldid=foldid,alpha=1,family=family,lower.limits=0,parallel=parallel)
  } else {
     meta.sta <- NULL
  }
  
  #--- simultaneous stacking ---
  if(any(stack=="sim")){
     meta.sim <- glmnet::cv.glmnet(y=y,x=cbind(y_hat[,1:k],X),alpha=alpha,family=family,
                               lower.limits=rep(c(0,-Inf),times=c(k,p)),
                               penalty.factor=rep(c(0,1),times=c(k,p)),foldid=foldid,
                               parallel=parallel)
  } else {
     meta.sim <- NULL
  }
  
  object <- list(base=base,meta.sta=meta.sta,meta.sim=meta.sim,prior.calib=prior.ext$beta,data=list(y=y,X=X,prior=prior),info=data.frame(n=n,p=p,k=k,family=family,alpha=alpha,scale=scale,stack=paste0(stack,collapse="+")))
  class(object) <- "transreg"
  
  return(object)
}

##### predict #####

#' @export
#' 
#' @title
#' Make Predictions
#' 
#' @description
#' Predicts outcome
#' 
#' @param object
#' object of class `transreg`
#' @inheritParams transreg
#' @param newx
#' features:
#' matrix with \eqn{n} rows (samples) and \eqn{p} columns (variables)
#' @param ... (not applicable)
#' 
#' @return
#' Returns predicted values or predicted probabilities.
#' The output is a column vector with one entry for each sample. 
#' 
#' @inherit transreg-package references
#' 
#' @inherit transreg seealso
#' 
#' @examples
#' #--- simulation ---
#' set.seed(1)
#' n0 <- 100; n1 <- 10000; n <- n0 + n1; p <- 500
#' X <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
#' beta <- rnorm(p)
#' prior <- beta + rnorm(p)
#' y <- X %*% beta
#' 
#' #--- train-test split ---
#' foldid <- rep(c(0,1),times=c(n0,n1))
#' y0 <- y[foldid==0]
#' X0 <- X[foldid==0,]
#' y1 <- y[foldid==1]
#' X1 <- X[foldid==1,]
#' 
#' #--- glmnet (without prior effects) ---
#' object <- glmnet::cv.glmnet(y=y0,x=X0)
#' y_hat <- predict(object,newx=X1,s="lambda.min")
#' mean((y1-y_hat)^2)
#' 
#' #--- transreg (with prior effects) ---
#' object <- transreg(y=y0,X=X0,prior=prior)
#' y_hat <- predict(object,newx=X1)
#' mean((y1-y_hat)^2) # decrease in MSE?
#' 
predict.transreg <- function(object,newx,stack=NULL,...){
  stack <- .which.stack(object,stack)
  eval(parse(text=paste0(".predict.",stack,"(object=object,newx=newx,...)")))
}

#' @title 
#' Internal functions
#'
#' @description
#' Internal functions called by
#' [coef.transreg()], [predict.transreg()] and [weights.transreg()], 
#' depending on choice between
#' standard stacking
#' and simultaneous stacking.
#'
#' @inheritParams predict.transreg
#'
#' @seealso
#' Use \code{\link[=coef.transreg]{coef}}, 
#' \code{\link[=predict.transreg]{predict}}
#' and \code{\link[=weights.transreg]{weights}}.
#'
#' @name extract
NULL

#' @describeIn extract called by `predict.transreg` if `stack="sta"`
.predict.sta <- function(object,newx,...){
  one <- newx %*% object$base$prior$beta
  two <- stats::predict(object$base,s=c(object$base$lambda.min,object$base$lambda.1se),newx=newx)
  y_hat <- stats::predict(object$meta.sta,s="lambda.min",newx=cbind(one,two),type="response")
  return(y_hat)
}

#' @describeIn extract called by `predict.transreg` if `stack="sim"`
.predict.sim <- function(object,newx,...){
  one <- newx %*% object$base$prior$beta
  y_hat <- stats::predict(object$meta.sim,s="lambda.min",newx=cbind(one,newx),type="response")
  return(y_hat)
}

##### coef #####

#' @export
#'
#' @title
#' Extract Coefficients
#'
#' @description
#' Extracts coefficients
#' from an object of class [transreg].
#'
#' @inheritParams predict.transreg
#'
#' @return
#' Returns estimated coefficients.
#' The output is a list with two slots:
#' slot `alpha` with the estimated intercept (scalar),
#' and slot `beta` with the estimated slopes (vector).
#'
#' @inherit transreg-package references
#' 
#' @inherit transreg seealso
#'
#' @examples
#' #--- simulation ---
#' set.seed(1)
#' n <- 100; p <- 500
#' X <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
#' beta <- rnorm(p)
#' prior <- beta + rnorm(p)
#' y <- X %*% beta
#' 
#' #--- glmnet (without prior effects) ---
#' object <- glmnet::cv.glmnet(y=y,x=X,alpha=0)
#' beta_hat <- coef(object,s="lambda.min")[-1]
#' mean((beta-beta_hat)^2)
#' 
#' #--- transreg (with prior effects) ---
#' object <- transreg(y=y,X=X,prior=prior,alpha=0)
#' beta_hat <- coef(object)$beta
#' mean((beta-beta_hat)^2) # decrease in MSE?
#'
coef.transreg <- function(object,stack=NULL,...){
  stack <- .which.stack(object,stack)
  eval(parse(text=paste0(".coef.",stack,"(object=object,...)")))
}

#' @describeIn extract called by `coef.transreg` if `stack="sta"`
.coef.sta <- function(object,...){
  beta <- stats::coef(object$base,s=c(object$meta.sta$lambda.min,object$meta.sta$lambda.1se))
  omega <- as.numeric(stats::coef(object$meta.sta,s=object$meta.sta$lambda.min))
  names <- paste0("source",seq_len(ncol(object$base$prior$beta)))
  names(omega) <- c("(Intercept)",names,"lambda.min","lambda.1se")
  p <- nrow(beta)-1
  alpha_star <- omega["(Intercept)"] + omega["lambda.min"]*beta["(Intercept)","s1"] + omega["lambda.1se"]*beta["(Intercept)","s2"]
  beta_star <- rep(NA,times=p)
  for(j in seq_len(p)){
    beta_star[j] <- sum(omega[names]*object$base$prior$beta[j,]) + omega["lambda.min"]*beta[1+j,"s1"] + omega["lambda.1se"]*beta[1+j,"s2"]
  }
  return(list(alpha=alpha_star,beta=beta_star))
}

#' @describeIn extract called by `coef.transreg` if `stack="sim"`
.coef.sim <- function(object,...){
  gamma <- object$base$prior$beta
  meta <- stats::coef(object$meta.sim,s="lambda.min")
  k <- ncol(object$base$prior$beta)
  p <- nrow(object$base$prior$beta)
  alpha_star <- meta[1]
  omega <- meta[2:(k+1)]
  betas <- meta[-c(1:(k+1))]
  beta_star <- rep(NA,times=p)
  for(i in seq_len(p)){
    beta_star[i] <- sum(omega*gamma[i,]) + betas[i]
  }
  return(list(alpha=alpha_star,beta=beta_star))
}

##### calibrate #####

#' @title 
#' Internal functions
#'
#' @description
#' Internal functions called by
#' [transreg()],
#' depending on choice between
#' exponential and isotonic calibration.
#'
#' @inheritParams transreg
#'
#' @seealso
#' Use [transreg()] for model fitting.
#'
#' @name calibrate
NULL

#' @describeIn calibrate called by `transreg` if `scale="exp"`
.exp.multiple <- function(y,X,prior,family,switch=FALSE,select=TRUE,track=FALSE){
  
  n <- nrow(X); p <- ncol(X); k <- ncol(prior)
  
  alpha <- theta <- pvalue <- sign <- tau <- rep(NA,times=k)
  beta <- matrix(NA,nrow=p,ncol=k)
  
  exp <- c(0,exp(seq(from=log(0.01),to=log(50),length.out=100)))
  
  for(j in seq_len(k)){
    coefs <- matrix(NA,nrow=length(exp),ncol=2)
    pred <- matrix(NA,nrow=nrow(X),ncol=length(exp))
    for(i in seq_along(exp)){
      temp <- sign(prior[,j])*abs(prior[,j])^exp[i]
      eta <- X %*% temp
      if(switch){
        glm <- stats::glm(y~eta,family=family)
        temp <- stats::coef(glm)
        coefs[i,1] <- temp["(Intercept)"]
        coefs[i,2] <- ifelse(is.na(temp["eta"]),0,temp["eta"])
        pred[,i] <- stats::fitted(glm)
      } else {
        if(stats::sd(eta)==0){eta <- eta + stats::rnorm(n,sd=1e-06); lambda <- 99e99}else{lambda <- 0}
        glmnet <- glmnet::glmnet(x=cbind(0,eta),y=y,family=family,lambda=lambda,lower.limits=0,intercept=TRUE)
        #temp <- stats::coef(glmnet)
        temp <- glmnet::coef.glmnet(glmnet)
        coefs[i,1] <- temp["(Intercept)","s0"]
        coefs[i,2] <- ifelse(is.na(temp["V2","s0"]),0,temp["V2","s0"])
        pred[,i] <- stats::predict(glmnet,newx=cbind(0,eta),s=lambda,type="response")
      }
    }
    cvm <- apply(pred,2,function(x) starnet:::.loss(y=y,x=x,family=family,type.measure="deviance"))
    id.min <- which.min(cvm)
    
    alpha[j] <- coefs[id.min,1]
    theta[j] <- coefs[id.min,2]
    beta[,j] <- coefs[id.min,2]*sign(prior[,j])*abs(prior[,j])^exp[id.min]
    sign[j] <- sign(coefs[id.min,2])
    tau[j] <- exp[id.min]
    
    if(select){
      # compare with empty model
      res0 <- .residuals(y=y,y_hat=mean(y),family=family)
      res1 <- .residuals(y=y,y_hat=pred[,which.min(cvm),drop=FALSE],family=family)
      pvalue[j] <- suppressWarnings(stats::wilcox.test(x=res1,y=res0,paired=TRUE,alternative="less")$p.value)
    }
    
    if(track){
      tryCatch(graphics::plot(x=exp,y=cvm,main=j),error=function(x) NULL)
      tryCatch(graphics::abline(v=exp[id.min]),error=function(x) NULL)
    }
  }
  
  if(select){
    remove <- pvalue > 0.05
    beta[,remove] <- 0
    if(track){message(ifelse(remove,".",ifelse(sign==1,"+","-")))}
  }
  
  return(list(alpha=alpha,beta=beta,theta=theta,tau=tau))
}

#' @describeIn calibrate called by `transreg` if `scale="iso"`
.iso.multiple <- function(y,X,prior,family,switch=FALSE,select=TRUE,track=FALSE){
  
  k <- ncol(prior)
  prior0 <- .iso.fast.single(y=y,X=X,prior=+prior,family=family)
  alpha0 <- prior0$alpha; beta0 <- prior0$beta
  fit0 <- joinet:::.mean.function(alpha0 + X %*% beta0,family=family)
  res <- .residuals(y=y,y_hat=mean(y),family=family)
  res0 <- .residuals(y=y,y_hat=fit0,family=family)
  pval0 <- apply(res0,2,function(x) suppressWarnings(stats::wilcox.test(x=x,y=res,paired=TRUE,alternative="less")$p.value))
  
  if(switch){
    prior1 <- .iso.fast.single(y=y,X=X,prior=-prior,family=family) # was .iso.slow.single
    alpha1 <- prior1$alpha; beta1 <- prior1$beta
    fit1 <- joinet:::.mean.function(alpha1 + X %*% beta1,family=family)
    res1 <- .residuals(y=y,y_hat=fit1,family=family)
    pval1 <- apply(res1,2,function(x) suppressWarnings(stats::wilcox.test(x=x,y=res,paired=TRUE,alternative="less")$p.value))
    if(track){message("p-value (+): ",paste0(signif(pval0,digits=2),sep=" "))}
    if(track){message("p-value (-): ",paste0(signif(pval1,digits=2),sep=" "))}
    cond <- pval0 <= pval1
    alpha0[!cond] <- alpha1[!cond]
    beta0[,!cond] <- beta1[,!cond]
    pvalue <- pmin(pval0,pval1)
  } else {
    cond <- rep(TRUE,times=k)
    pvalue <- pval0
  }
  
  if(select){
    remove <- pvalue > 0.05
    beta0[,remove] <- 0
    if(track){message(ifelse(remove,".",ifelse(cond,"+","-")))}
  }
  
  return(list(alpha=alpha0,beta=beta0))
}

#' @describeIn calibrate called by `transreg` if `scale="iso"` (via `.iso.multiple`)
.iso.fast.single <- function(y,X,prior,family){
  n <- length(y)
  p <- nrow(prior)
  k <- ncol(prior)
  ALPHA <- rep(NA,times=k)
  BETA <- posterior <- matrix(data=NA,nrow=p,ncol=k)
  for(i in seq_len(k)){
    single.prior <- prior[,i]
    order <- order(single.prior)
    Xo <- X[,order,drop=FALSE]
    sign <- sign(single.prior[order])
    A <- t(apply(Xo[,sign<0,drop=FALSE],1,function(x) cumsum(x)))
    B <- t(apply(Xo[,sign>=0,drop=FALSE],1,function(x)  rev(cumsum(rev(x)))))
    if(any(sign<0)&any(sign>=0)){
      if(nrow(A)==1){A <- matrix(A,ncol=1)}
      if(nrow(B)==1){B <- matrix(B,ncol=1)}
      Xcs <- cbind(A,B)
    } else {
      if(any(sign<0)){
        Xcs <- A
      } else if(any(sign>=0)){
        Xcs <- B
      }
    }
    times <- c(sum(sign==-1),sum(sign==0),sum(sign==+1))
    inc <- glmnet::glmnet(y=y,x=Xcs,family=family,alpha=0,lambda=0,lower.limits=rep(c(-Inf,0,0),times=times),upper.limits=rep(c(0,0,Inf),times=times))
    coef <- stats::coef(inc)
    ALPHA[i] <- coef[1]
    beta <- coef[-1]
    new <- c(rev(cumsum(rev(beta[sign<0]))),cumsum(beta[sign>=0]))
    BETA[,i] <- new[order(order)]
  }
  return(list(alpha=ALPHA,beta=BETA))
}

#' @describeIn calibrate replaced by `.iso.fast.single`
.iso.slow.single <- function(y,X,prior,family){
  
  n <- length(y)
  p <- nrow(prior)
  k <- ncol(prior)
  
  prior[-1e-06 < prior & prior < +1e-06] <- 0
  
  ALPHA <- rep(NA,times=k)
  BETA <- posterior <- matrix(data=NA,nrow=p,ncol=k)
  for(i in seq_len(k)){
    single.prior <- prior[,i]
    p <- ncol(X)
    order <- order(single.prior)
    prior.order <- single.prior[order]
    X.order <- X[,order]
    alpha <- CVXR::Variable(1)
    beta <- CVXR::Variable(p)
    constraints <- list()
    constraints$iso <- CVXR::diff(beta)>=0
    if(any(prior.order<0)){
      constraints$neg <- beta[prior.order<0]<=0
    }
    if(any(prior.order==0)){
      constraints$zero <- beta[prior.order==0]==0
    }
    if(any(prior.order>0)){
      constraints$pos <- beta[prior.order>0]>=0
    }
    if(family=="gaussian"){
      loss <- CVXR::sum_squares(y - alpha - X.order %*% beta) / (2 * n)
    } else if(family=="binomial"){
      loss <- sum(alpha+X.order[y==0,]%*%beta)+sum(CVXR::logistic(-alpha-X.order%*%beta))
    } else {
      stop("Family not yet implemented!")
    }
    penalty <- (1e-06)*CVXR::p_norm(beta,p=2) # lasso: p=1, ridge: p=2 # was p=2 (ridge)
    
    # The ridge penalty should be multiplied by a very small regularisation parameter (e.g. 1e-6). Why is this insufficient for the binomial case? (For the binomial case, an alternative might be to impose a lower limit of minus one and an upper limit of plus one on the coefficients.)
    
    objective <- CVXR::Minimize(loss+penalty)
    problem <- CVXR::Problem(objective=objective,constraints=constraints)
    result <- CVXR::solve(problem)
    
    if(result$status=="solver_error"){
      warning("Solver error ...")
      ALPHA[i] <- mean(y)
      BETA[,i] <- 0
    } else {
      ALPHA[i] <- result$getValue(alpha)
      beta.r <- result$getValue(beta)
      BETA[,i] <- beta.r[order(order)]
    }
  }
  
  BETA[-1e-06 < BETA & BETA < +1e-06] <- 0
  return(list(alpha=ALPHA,beta=BETA))
}

##### reproducibility #####

#' @title 
#' Cross-validation (reproducibility)
#' 
#' @description
#' Function for reproducing hold-out method (simulation)
#' and \eqn{k}-fold cross-validation (application).
#' See vignette.
#' 
#' @inheritParams transreg
#' @param target
#' list with slot x (feature matrix with n rows and p columns) and slot y (target vector of length n)
#' @param source
#' list of k lists, each with slot x (feature matrix with m_i rows and p columns) and slot y (target vector of length m_i)
#' @param alpha.prior
#' number between 0 (lasso) and 1 (ridge), character "p-value", or NULL (alpha.prior=alpha, but if alpha=1 then alpha.prior=0.95)
#' @param z
#' prior weights
#' @param foldid.ext
#' external fold identifiers
#' @param nfolds.ext
#' number of external folds
#' @param foldid.int
#' internal fold identifiers
#' @param nfolds.int
#' number of internal folds
#' @param type.measure
#' character
#' @param alpha.prior
#' alpha for source regression
#' @param naive
#' compare with naive transfer learning: logical
#' @param seed
#' random seed
#' @param cores
#' number of cores for parallel computing
#' (requires R package `doMC`)
#' @param xrnet
#' compare with xrnet: logical
#' 
#' @seealso 
#' [transreg()]
compare <- function(target,source=NULL,prior=NULL,z=NULL,family,alpha,scale="iso",sign=FALSE,switch=FALSE,select=TRUE,foldid.ext=NULL,nfolds.ext=10,foldid.int=NULL,nfolds.int=10,type.measure="deviance",alpha.prior=NULL,naive=TRUE,seed=NULL,cores=1,xrnet=FALSE){
  
  if(cores>1){
    doMC::registerDoMC(cores=cores)
  }
  
  if(!is.null(z) && any(z<0)){
    stop("Prior weights (argument \"z\") must be non-negative.")
  }
  
  if(is.null(source)==is.null(prior)){
    stop("Provide either \"source\" or \"prior\".",call.=FALSE)
  }
  
  if(!is.list(target)||is.null(names(target))){stop("Expect argument target as list with slots x and y.")}
  if(!any(names(target) %in% c("y","x"))){stop("Expect argument target as list with slots x and y.")}
  p <- ncol(target$x)
  
  if(!is.null(source)){
  if(!is.list(source)){stop("Expect argument source as list of lists.")}
  for(i in seq_along(source)){
    if(!is.list(source[[i]])||is.null(names(source[[i]]))){stop("Expect argument source as list of lists with slots x and y.")}
    if(!any(names(source[[i]]) %in% c("y","x"))){stop("Expect argument source as list of lists with slots x and y.")}
  }
  }
  
  if(is.vector(prior)){
    prior <- matrix(prior,ncol=1)
  }
  
  if(is.data.frame(target$x)){
    target$x <- as.matrix(target$x)
  }

  if(!is.null(source)){
    target$x <- scale(target$x)
    target$x[is.na(target$x)] <- 0
    for(i in seq_along(source)){
      if(is.data.frame(source[[i]]$x)){
        source[[i]]$x <- as.matrix(source[[i]]$x)
      }
      source[[i]]$x <- scale(source[[i]]$x)
      source[[i]]$x[is.na(source[[i]]$x)] <- 0
    }
    
      # # alternative standardisation
      # temp <- rbind(target$x,do.call(what="rbind",args=lapply(source,function(x) x$x)))
      # mu_temp <- colMeans(temp)
      # sd_temp <- apply(temp,2,sd)
      # n <- c(length(target$y),sapply(source,function(x) length(x$y)))
      # mus <- rbind(colMeans(target$x),t(sapply(source,function(x) colMeans(x$x))))
      # mu <- colSums(mus*n)/sum(n)
      # sd <- sqrt(rowSums(cbind(rowSums((t(target$x)-mu)^2),sapply(source,function(x) rowSums((t(x$x)-mu)^2))))/(sum(n)-1))
      # if(any(abs(mu-mu_temp)>1e-06)){stop("Invalid.")}
      # if(any(abs(sd-sd_temp)>1e-06)){stop("Invalid.")}
      # 
      # target$x <- t((t(target$x)-mu)/sd)
      # target$x[,sd==0] <- 0
      # for(i in seq_along(source)){
      #   source[[i]]$x <- t((t(source[[i]]$x)-mu)/sd)
      #   source[[i]]$x[,sd==0] <- 0
      # }
    
    if(is.null(alpha.prior)){alpha.prior <- ifelse(alpha==1,0.95,alpha)}
    
    # prior
    k <- length(source)
    prior <- matrix(NA,nrow=p,ncol=k)
    for(i in seq_len(k)){
      foldid.source <- .folds(y=source[[i]]$y,nfolds=10)$foldid
      if(is.character(alpha.prior) & alpha.prior=="p-value"){
        #if(family=="binomial"){
        #  p.value <- apply(source[[i]]$x,2,function(x) suppressWarnings(stats::wilcox.test(x~source[[i]]$y)$p.value))
        #  sign <- apply(source[[i]]$x,2,function(x) mean(x[source[[i]]$y==1])-mean(x[source[[i]]$y==0]))
        #  prior[,i] <- sign*(-log10(p.value))
        #} else {
        #  stop("Implement correlation test!",call.=FALSE)
        #}
        for(j in seq_len(p)){
          yy <- source[[i]]$y
          xx <- source[[i]]$x[,j]
          lm <- summary(stats::glm(yy~xx,family=family))$coefficients
          if(nrow(lm)==1 | any(is.na(lm))){prior[j,i] <- 0; next}
          sign <- sign(lm["xx","Estimate"])
          p.value <- lm["xx",ifelse(family=="gaussian","Pr(>|t|)","Pr(>|z|)")]
          prior[j,i] <- sign*(-log10(p.value))
        }
      } else {
        object <- glmnet::cv.glmnet(y=source[[i]]$y,x=source[[i]]$x,family=family,alpha=alpha.prior,foldid=foldid.source,parallel=cores>1)
        prior[,i] <- as.numeric(stats::coef(object,s="lambda.min"))[-1]
      }
    }
  } else {
    k <- ncol(prior)
  }
  
  if(!is.null(seed)){set.seed(seed)}
  folds <- .folds(y=target$y,nfolds.ext=nfolds.ext,nfolds.int=nfolds.int,foldid.ext=foldid.ext,foldid.int=foldid.int)
  foldid.ext <- folds$foldid.ext
  foldid.int <- folds$foldid.int
  
  temp <- paste0(rep(c("transreg.","transreg."),each=length(scale)),scale,rep(c(".sta",".sim"),each=length(scale)))
  names <- c("mean","glmnet","glmtrans"[!is.null(source)],temp,"fwelnet"[!is.null(z)],"ecpc"[!is.null(z)],"naive"[naive],"xrnet"[xrnet])
  pred <- matrix(data=NA,nrow=length(target$y),ncol=length(names),dimnames=list(NULL,names))
  coef <- sapply(names,function(x) matrix(data=NA,nrow=ncol(target$x),ncol=nfolds.ext),simplify=FALSE)
  time <- rep(0,time=length(names)); names(time) <- names
  time <- time[!grepl(pattern="transreg",x=names(time))]
  time["transreg"] <- 0
  
  for(i in seq_len(nfolds.ext)){
    y0 <- target$y[foldid.ext!=i]
    X0 <- target$x[foldid.ext!=i,]
    X1 <- target$x[foldid.ext==i,]
    
    foldid <- foldid.int[foldid.ext!=i]
    
    # mean
    if(!is.null(seed)){set.seed(seed)}
    start <- Sys.time()
    pred[foldid.ext==i,"mean"] <- rep(mean(y0),times=sum(foldid.ext==i))
    coef$mean[,i] <- rep(0,times=p)
    end <- Sys.time()
    time["mean"] <- time["mean"]+difftime(time1=end,time2=end,units="secs")
    
    # glmnet
    if(!is.null(seed)){set.seed(seed)}
    start <- Sys.time()
    object <- glmnet::cv.glmnet(y=y0,x=X0,family=family,foldid=foldid,alpha=alpha,parallel=cores>1)
    end <- Sys.time()
    pred[foldid.ext==i,"glmnet"] <- as.numeric(stats::predict(object,newx=X1,s="lambda.min",type="response"))
    coef$glmnet[,i] <- coef(object,s="lambda.min")[-1]
    time["glmnet"] <-  time["glmnet"]+difftime(time1=end,time2=start,units="secs")
    
    # glmtrans
    if(!is.null(source)){
      if(!is.null(seed)){set.seed(seed)}
      start <- Sys.time()
      object <- tryCatch(glmtrans::glmtrans(target=list(x=X0,y=y0),source=source,alpha=alpha,family=family,nfolds=nfolds.int),error=function(x) NULL)
      end <- Sys.time()
      if(!is.null(object)){
        pred[foldid.ext==i,"glmtrans"] <- stats::predict(object,newx=X1,type="response")
        coef$glmtrans[,i] <- object$beta[-1]
      }
      time["glmtrans"] <-  time["glmtrans"]+difftime(time1=end,time2=start,units="secs")
    }
    
    # transreg
    for(j in seq_along(scale)){
      if(!is.null(seed)){set.seed(seed)}
      start <- Sys.time()
      object <- transreg(y=y0,X=X0,prior=prior,family=family,foldid=foldid,alpha=alpha,scale=scale[j],stack=c("sta","sim"),sign=sign,switch=switch,select=select,parallel=cores>1)
      end <- Sys.time()
      pred[foldid.ext==i,paste0("transreg.",scale[j],".sta")] <- .predict.sta(object,newx=X1)
      pred[foldid.ext==i,paste0("transreg.",scale[j],".sim")] <- .predict.sim(object,newx=X1)
      coef[[paste0("transreg.",scale[j],".sta")]][,i] <- coef(object,stack="sta")$beta
      coef[[paste0("transreg.",scale[j],".sim")]][,i] <- coef(object,stack="sim")$beta
      time["transreg"] <- time["transreg"]+difftime(time1=end,time2=start,units="secs")
    }
    
    # naive
    if(naive){
      if(!is.null(seed)){set.seed(seed)}
      start <- Sys.time()
      glm <- stats::glm(formula=y0 ~ .,family=family,data=data.frame(x=X0 %*% prior)) # was y~x
      end <- Sys.time()
      pred[foldid.ext==i,"naive"] <- stats::predict(glm,newdata=data.frame(x=X1 %*% prior),type="response")
      time["naive"] <- time["naive"] + difftime(time1=end,time2=start,units="secs")
    }
    
    if(xrnet){
      if(!is.null(seed)){set.seed(seed)}
      start <- Sys.time()
      cond <- apply(prior,2,stats::sd)>0
      tune_xrnet <- xrnet::tune_xrnet(x=X0,y=y0,
                                      external=prior[,cond,drop=FALSE],
                                      penalty_main=xrnet::define_penalty(penalty_type=alpha),
                                      family=family,
                                      foldid=foldid,
                                      nfolds=nfolds.int,
                                      loss="deviance")
      end <- Sys.time()
      pred[foldid.ext==i,"xrnet"] <- stats::predict(tune_xrnet,newdata=X1) # was xrnet:::predict.tune_xrnet
      time["xrnet"] <- time["xrnet"] + difftime(time1=end,time2=start,units="secs")
    }
    
    # co-data methods
    if(!is.null(z)){
      
      # fwelnet
      if(!is.null(seed)){set.seed(seed)}
      start <- Sys.time()
      object <- tryCatch(fwelnet::cv.fwelnet(x=X0,y=y0,z=z,family=family,foldid=foldid,alpha=alpha),error=function(x) NULL)
      end <- Sys.time()
      if(!is.null(object)){
        pred[foldid.ext==i,"fwelnet"] <- stats::predict(object,xnew=X1,s="lambda.min",type="response")
        coef$fwelnet[,i] <- object$glmfit$beta[,which.min(object$cvm)]
      }
      time["fwelnet"] <- time["fwelnet"]+difftime(time1=end,time2=start,units="secs")

      # ecpc - default
      if(!is.null(seed)){set.seed(seed)}
      start <- Sys.time()
      if(is.vector(z)){
        Z <- list(as.numeric(z))
      } else {
        Z <- as.list(as.data.frame(z))
      }
      #--- original ---
      #object <- tryCatch(ecpc::ecpc(Y=y0,X=X0,Z=Z,X2=X1,fold=nfolds.int),error=function(x) NULL)
      #--- splines start ---
      co.Z <- co.S <- co.C <- list()
      for(j in seq_len(k)){
        name <- paste0("Z",j)
        co.Z[[name]] <- ecpc::createZforSplines(values=Z[[j]])
        co.S[[name]] <- list(ecpc::createS(orderPen=2,G=ncol(co.Z[[name]])))
        co.C[[name]] <- ecpc::createCon(G=ncol(co.Z[[name]]),shape="positive+monotone.i")
      }
      object <- tryCatch(ecpc::ecpc(Y=y0,X=X0,Z=co.Z,paraPen=co.S,paraCon=co.C,X2=X1,fold=nfolds.int),error=function(x) NULL)
      end <- Sys.time()
      #--- splines end ---
      if(!is.null(object)){
        pred[foldid.ext==i,"ecpc"] <- object$Ypred
        coef$ecpc[,i] <- coef(object)$beta
      }
      time["ecpc"] <- time["ecpc"]+difftime(time1=end,time2=start,units="secs")
      
      # #--- GRridge ---
      # if(is.null(partitions)|is.null(monotone)){stop()}
      # if(!is.null(seed)){set.seed(seed)}
      # start <- Sys.time()
      # object <- tryCatch(GRridge::grridge(highdimdata=t(X0),response=y0,
      #                                       partitions=partitions,monotone=monotone,
      #                                       innfold=nfolds.int,fixedfoldsinn=TRUE),
      #                      error=function(x) NULL)
      # if(!is.null(object)){
      #   temp <- GRridge::predict.grridge(object=object,datanew=as.data.frame(t(X1)))
      #   pred[foldid.ext==i,"NoGroups"] <- temp[,"NoGroups"]
      #   pred[foldid.ext==i,"GRridge"] <- temp[,"GroupRegul"]
      # }
      # end <- Sys.time()
      # time["GRridge"] <- time["GRridge"]+difftime(time1=end,time2=start,units="secs")
      
      # #--- xtune ---
      # if(!is.null(seed)){set.seed(seed)}
      # start <- Sys.time()
      # model <- switch(family,"gaussian"="linear","binomial"="binary")
      # method <- ifelse(alpha==0,"ridge",ifelse(alpha==1,"lasso",stop()))
      # cond <- !duplicated(t(X0))
      # object <- tryCatch(xtune::xtune(Y=y0,X=X0[,cond],Z=z[cond,],family=model,method=method),error=function(x) NULL)
      # if(!is.null(object)){
      #   temp <- stats::predict(object,newX=X1[,cond],type="response")
      #   #temp <- joinet:::.mean.function(temp,family=family) # trial!
      #   # It is possible that xtune 0.10 returns the linear predictors and not the predicted probabilities!
      #   pred[foldid.ext==i,"xtune"] <- temp
      # }
      # end <- Sys.time()
      # time["xtune"] <- time["xtune"]+difftime(time1=end,time2=start,units="secs")

      # #--- CoRF ---
      # if(!is.null(seed)){set.seed(seed)}
      # CoData <- as.data.frame(z)
      # CoDataRelation <- rep("increasing",times=ncol(CoData))
      # start <- Sys.time()
      # object <- tryCatch(CoRF::CoRF(Y=y0,X=X0,CoData=CoData,
      #               CoDataRelation=CoDataRelation),
      #               error=function(x) NULL)
      # if(!is.null(object)){
      #   temp <- stats::predict(object$SavedForests[[length(object$SavedForests)]],newdata=data.frame(Xdf=X1))
      #   pred[foldid.ext==i,"CoRF"] <- temp$predicted[,2]
      # }
      # end <- Sys.time()
      # time["CoRF"] <- time["CoRF"]+difftime(time1=end,time2=start,units="secs")

    }
    
  }
  
  pred <- pred[,colMeans(is.na(pred))!=1]
  loss <- list()
  for(i in seq_along(type.measure)){
    loss[[type.measure[i]]] <- apply(pred,2,function(x) starnet:::.loss(y=target$y[foldid.ext!=0],x=x[foldid.ext!=0],
                                              family=family,type.measure=type.measure[i],foldid=foldid.int[foldid.ext!=0]))
  }
  
  list <- loss
  list$time <- time
  list$coef <- coef
  list$foldid <- foldid.ext
  list$pred <- pred
  list$p <- p
  
  return(list)
}

# @example (if function 'compare' gets exported)
# set.seed(1)
# n <- 100; p <- 500
# X <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
# beta <- rnorm(p)*rbinom(n=p,size=1,prob=0.2)
# y <- X %*% beta
# \dontshow{
# subset <- 1:10 # speed up for CRAN
# object <- suppressMessages(transreg:::compare(target=list(y=y,x=X[,subset]),prior=beta[subset],family="gaussian",alpha=0))}
# \donttest{
# object <- transreg:::compare(target=list(y=y,x=X),prior=beta,family="gaussian",alpha=0)}

#' @title
#' Simulation (reproducibility)
#' 
#' @description 
#' Function for reproducing 'internal' simulation study.
#' See vignette.
#' 
#' @param p
#' number of features
#' @param n.target
#' sample size for target data set
#' @param n.source
#' sample size(s) for source data set(s), scalar or vector of length k
#' @param family
#' "Gaussian", "binomial" or "poisson"
#' @param k
#' number of source data sets
#' @param prop
#' approximate proportion of features with effects
#' @param rho.beta
#' correlation between effects (across different data sets)
#' @param rho.x
#' base for decreasing correlation structure for correlation between features
#' @param w
#' weight between signal and noise
#' @param trans
#' logical vector of length \eqn{k}:
#' transferable (TRUE) or non-transferable (FALSE) source
#' @param exp
#' non-negative vector of length \eqn{k}
#' for transforming beta to sign(beta)*abs(beta)^exp
#' 
#' @seealso
#' Use [glmtrans::models()] for reproducing 'external' simulation study.
#' 
simulate <- function(p=1000,n.target=100,n.source=150,k=2,family="gaussian",prop=0.01,rho.beta=0.95,rho.x=0.95,w=0.5,trans=rep(TRUE,times=k),exp=rep(1,times=k)){
  
  target <- source <- list()
  
  # effects
  mu <- rep(x=0,times=k+1)
  Sigma <- matrix(data=rho.beta,nrow=k+1,ncol=k+1)
  Sigma[c(!trans,FALSE),] <- 0; Sigma[,c(!trans,FALSE)] <- 0
  diag(Sigma) <- 1
  beta <- mvtnorm::rmvnorm(n=p,mean=mu,sigma=Sigma)
  
  #message("This part is temporary!")
  for(i in seq_len(k)){
    beta[,i] <- sign(beta[,i])*abs(beta[,i])^exp[i]
  }
  
  cond <- mvtnorm::rmvnorm(n=p,mean=mu,sigma=Sigma)>stats::qnorm(1-prop)
  beta[cond==0] <- 0
  
  # features
  mu <- rep(x=0,times=p)
  Sigma <- matrix(data=NA,nrow=p,ncol=p)
  Sigma <- rho.x^abs(row(Sigma)-col(Sigma))
  diag(Sigma) <- 1
  
  # source
  for(i in seq_len(k)){
    source[[i]] <- list()
    source[[i]]$x <- mvtnorm::rmvnorm(n=n.source,mean=mu,sigma=Sigma)
    eta <- sqrt(w)*as.vector(scale(source[[i]]$x %*% beta[,i])) + sqrt(1-w)*stats::rnorm(n.source)
    source[[i]]$y <- joinet:::.mean.function(eta,family=family)
    if(family %in% c("binomial","poisson")){source[[i]]$y <- round(source[[i]]$y)}
  }
  
  # target
  target$x <- mvtnorm::rmvnorm(n=n.target,mean=mu,sigma=Sigma)
  eta <- sqrt(w)*as.vector(scale(target$x %*% beta[,k+1])) + sqrt(1-w)*stats::rnorm(n.target)
  target$y <- joinet:::.mean.function(eta,family=family)
  if(family %in% c("binomial","poisson")){target$y <- round(target$y)}
  
  return(list(source=source,target=target,beta=beta))
}

## example (if function 'simulate' gets exported)
# set.seed(1)
# data <- transreg:::simulate(p=200)
# prior <- numeric()
# for(i in seq_along(data$source)){
#   glmnet <- glmnet::cv.glmnet(y=data$source[[i]]$y,x=data$source[[i]]$x)
#   prior <- cbind(prior,coef(glmnet,s="lambda.min")[-1])
#}
# \donttest{
# object <- transreg(y=data$target$y,X=data$target$x,prior=prior)}
# 

##### weights #####

#' @export
#' @importFrom stats weights
#' 
#' @title
#' Extract Weights
#'
#' @description
#' Extracts weights from an object of class [transreg].
#' 
#' @inheritParams predict.transreg
#'
#' @return
#' Returns weights.
#' The output is a numerical vector
#' with one entry for each source of co-data.
#'
#' @inherit transreg-package references
#' 
#' @seealso
#' This function is about weights for sources of prior effects.
#' To extract weights for features (estimated regression coefficients),
#' use [coef()].
#'
#' @examples
#' #--- simulation ---
#' set.seed(1)
#' n <- 100; p <- 500
#' X <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
#' beta <- rnorm(p)
#' prior <- cbind(beta+rnorm(p),beta+rnorm(p),rnorm(p),rnorm(p))
#' y <- X %*% beta
#' 
#' #--- example ---
#' object <- transreg(y=y,X=X,prior=prior)
#' weights(object)
#'
weights.transreg <- function(object,stack=NULL,...){
  stack <- .which.stack(object,stack)
  eval(parse(text=paste0(".weights.",stack,"(object=object,...)")))
}

#' @describeIn extract called by `weights.transreg` if `stack="sta"`
.weights.sta <- function(object,...){
  stats::coef(object$meta.sta,s="lambda.min")[2:(object$info$k+1)]
}

#' @describeIn extract called by `weights.transreg` if `stack="sim"`
.weights.sim <- function(object,...){
  stats::coef(object$meta.sim,s="lambda.min")[2:(object$info$k+1)]
}

##### other #####

#' @export
#' 
#' @title 
#' Print transreg-object
#' 
#' @description
#' Show summary of transreg-object
#' 
#' @param x object of class transreg
#' @inheritParams predict.transreg
#' 
#' @return
#' Returns family of distributions,
#' elastic net mixing parameter (\eqn{alpha}),
#' number of samples (\eqn{n}),
#' number of features (\eqn{p}),
#' number of sources of co-data (\eqn{k}),
#' chosen calibration method (exponential or isotonic),
#' and chosen stacking method (standard or simultaneous).
#' 
#' @examples
#' #--- simulation ---
#' set.seed(1)
#' n <- 100; p <- 500
#' X <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
#' beta <- rnorm(p)
#' prior <- beta + rnorm(p)
#' y <- X %*% beta
#' 
#' #--- print.transreg  ---
#' object <- transreg(y=y,X=X,prior=prior)
#' object
#' 
print.transreg <- function(x,...){
  cat("----- transreg-object -----\n")
  cat(paste0("family: '",x$info$family,"'\n"))
  name <- ifelse(x$info$alpha==0,"ridge",ifelse(x$info$alpha==1,"lasso","elastic net"))
  cat(paste0("alpha = ",x$info$alpha," (",name,")\n"))
  cat("n =",x$info$n,"(samples)\n")
  cat("p =",x$info$p,"(features)\n")
  cat("k =",x$info$k,"(sources of co-data)\n")
  cat(paste0("calibration: '",x$info$scale,"'\n"))
  names <- c("sta"[!is.null(x$meta.sta)],"sim"[!is.null(x$meta.sim)])
  cat("stacking:",paste(paste0("'",names,"'"),collapse=" and "),"\n")
  cat("---------------------------")
}

#' @export
#' @importFrom stats fitted
#' 
#' @title
#' Fitted values
#'
#' @description
#' Extracts fitted values
#' 
#' @inheritParams predict.transreg
#'
#' @return
#' Returns fitted values.
#' The output is a numerical vector
#' with one entry for sample.
#'
#' @inherit transreg-package references
#' 
#' @inherit transreg seealso
#'
#' @examples
#' #--- simulation ---
#' set.seed(1)
#' n0 <- 100; n1 <- 10000; n <- n0 + n1; p <- 500
#' X <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
#' beta <- rnorm(p)
#' prior <- beta + rnorm(p)
#' y <- X %*% beta
#' 
#' #--- train-test split ---
#' foldid <- rep(c(0,1),times=c(n0,n1))
#' y0 <- y[foldid==0]
#' X0 <- X[foldid==0,]
#' y1 <- y[foldid==1]
#' X1 <- X[foldid==1,]
#' 
#' object <- transreg(y=y0,X=X0,prior=prior)
#' 
#' #--- fitted values ---
#' y0_hat <- fitted(object)
#' mean((y0-y0_hat)^2)
#' 
#' #--- predicted values ---
#' y1_hat <- predict(object,newx=X1)
#' mean((y1-y1_hat)^2) # increase in MSE?
#' 
fitted.transreg <- function(object,stack=NULL,...){
  stats::predict(object,newx=object$data$X,stack=stack,...)
}

#' @export
#' 
#' @title
#' Plot transreg-object
#'
#' @description
#' Plot transreg-object
#' 
#' @param x object of type transreg
#' @inheritParams predict.transreg
#'
#' @return
#' Returns four plots.
#' 
#' * top-left:
#' Calibrated prior effects (\eqn{y}-axis) against
#' original prior effects (\eqn{x}-axis).
#' Each line is for one source of prior effects,
#' with the colour given by [grDevices::palette()]
#' (black: 1, red: 2, green: 3, blue: 4, ...).
#' 
#' * top-right:
#' Estimated coefficients with transfer learning (\eqn{y}-axis)
#' against estimated coefficients without transfer learning (\eqn{x}-axis).
#' Each point represents one feature.
#' 
#' * bottom-left:
#' Estimated weights for sources of prior effects
#' (labels 1 to \eqn{k}),
#' and either
#' estimated weights for `lambda.min` and `lambda.1se` models
#' (standard stacking)
#' or estimated weights for features
#' (simultaneous stacking).
#' 
#' * bottom-right:
#' Absolute deviance residuals (\eqn{y}-axis)
#' against fitted values (\eqn{x}-axis).
#' Each point represents one sample.
#'
#' @inherit transreg-package references
#' 
#' @inherit transreg seealso
#'
#' @examples
#' #--- simulation ---
#' set.seed(1)
#' n <- 100; p <- 500
#' X <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
#' beta <- rnorm(p) #*rbinom(n=n,size=1,prob=0.2)
#' prior1 <- beta + rnorm(p)
#' prior2 <- beta + rnorm(p)
#' prior3 <- rnorm(p)
#' prior4 <- rnorm(p)
#' y <- X %*% beta
#' 
#' prior <- cbind(prior1,prior2,prior3,prior4)
#' object <- transreg(y=y,X=X,prior=prior,alpha=0,stack=c("sta","sim"))
#' 
#' plot(object,stack="sta")
#' 
plot.transreg <- function(x,stack=NULL,...){
  object <- x
  stack <- .which.stack(object,stack=stack)
  scale <- switch(object$info$scale,"exp"="exponential","iso"="isotonic","?")
  oldpar <- graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow=c(2,2),mar=c(3,3,1,1))

  #--- calibrated vs initial prior effects ---
  graphics::plot.new()
  graphics::plot.window(xlim=range(object$data$prior),ylim=range(object$prior.calib))
  graphics::box()
  graphics::axis(side=1,cex.axis=0.8)
  graphics::axis(side=2,cex.axis=0.8)
  graphics::title(main=paste(scale,"calibration"),cex.main=0.8,line=0.2)
  graphics::title(xlab="initial prior",ylab="calibrated prior",line=2,cex.lab=0.8)
  graphics::abline(h=0,col="grey",lty=2)
  graphics::abline(v=0,col="grey",lty=2)
  for(i in seq_len(object$info$k)){
    x <- object$data$prior[,i]
    y <- object$prior.calib[,i]
    graphics::lines(x=x[order(x)],y=y[order(x)],col=i)
  }
  
  #--- estimated betas without and with prior effects ---
  x <- stats::coef(object$base,s="lambda.min")[-1]
  y <- stats::coef(object,stack=stack)$beta
  graphics::plot.new()
  graphics::plot.window(xlim=range(x),ylim=range(y))
  graphics::box()
  graphics::axis(side=1,cex.axis=0.8)
  graphics::axis(side=2,cex.axis=0.8)
  graphics::title(main="estimated effects",cex.main=0.8,line=0.2)
  graphics::title(xlab="without prior",ylab="with prior",line=2,cex.lab=0.8)
  graphics::abline(h=0,col="grey",lty=2)
  graphics::abline(v=0,col="grey",lty=2)
  graphics::points(x=x,y=y,cex=0.5)
  lm <- lm(y~x)
  graphics::lines(x=x,y=fitted(lm),col="grey")
  graphics::legend(x="bottomright",legend=paste0("p=",object$info$p),bty="n",cex=0.8)
  
  #--- weights for sources of prior effects ---
  x <- seq_len(object$info$k)
  y <- weights(object,stack=stack)
  if(stack=="sta"){
    z <- stats::coef(object$meta.sta,s="lambda.min")[-c(1:(object$info$k+1))]
  } else if(stack=="sim"){
    z <- abs(stats::coef(object$meta.sim,s="lambda.min")[-c(1:(object$info$k+1))])
  } else {
    stop("Invalid.")
  }
  graphics::plot.new()
  graphics::plot.window(xlim=c(0.5,object$info$k+2),ylim=c(0,max(c(y,z))))
  graphics::box()
  graphics::axis(side=1,cex.axis=0.8,at=seq_len(object$info$k))
  graphics::axis(side=2,cex.axis=0.8)
  graphics::title(main="sources",cex.main=0.8,line=0.2)
  graphics::title(xlab="source",ylab="weight",line=2,cex.lab=0.8)
  graphics::segments(x0=x,y0=0,y1=y,lwd=3,col=x)
  graphics::abline(v=object$info$k+0.5,col="grey",lty=3)
  graphics::segments(x0=seq(from=object$info$k+1,to=object$info$k+2,length.out=length(z)),y0=0,y1=z,lwd=ifelse(stack=="sta",3,1))
  if(stack=="sta"){
    graphics::axis(side=1,cex.axis=0.8,at=object$info$k+c(1,2),label=c("min","1se"))
  } else if(stack=="sim"){
    graphics::axis(side=1,cex.axis=0.8,at=object$info$k+1.5,label="features",tick=FALSE)
  }
  
  #--- outliers ---
  y_hat <- stats::fitted(object,stack=stack)
  resid <- .residuals(y=object$data$y,y_hat=y_hat,family=object$info$family)
  graphics::plot.new()
  graphics::plot.window(xlim=range(y_hat),ylim=range(resid))
  graphics::box()
  graphics::axis(side=1,cex.axis=0.8)
  graphics::axis(side=2,cex.axis=0.8)
  graphics::title(main="samples",cex.main=0.8,line=0.2)
  xlab <- paste("fitted",ifelse(object$info$family=="binomial","probability","value"))
  graphics::title(xlab=xlab,ylab="deviance",line=2,cex.lab=0.8)
  graphics::abline(h=0,col="grey",lty=2)
  graphics::points(x=y_hat,y=resid,cex=0.5)
  graphics::legend(x="bottomright",legend=paste0("n=",object$info$n),bty="n",cex=0.8)

}

#' @title Calculate residuals
#' 
#' @description
#' Calculates residuals from observed outcome
#' and predicted values (Gaussian family)
#' or predicted probabilities (binomial family).
#' Called by `.exp.multiple` and `.iso.multiple`.
#' 
#' @param y
#' response: vector of length \eqn{n} (see family)
#' @param y_hat
#' predicted values or probabilities (see family):
#' vector of length \eqn{n},
#' or matrix with \eqn{n} rows (samples) and \eqn{k} columns (methods)
#' @param family
#' character
#' "gaussian" (\eqn{y}: real numbers, \eqn{y\_hat}: real numbers)
#' or "binomial" (\eqn{y}: 0s and 1s, \eqn{y\_hat}: unit interval)
#' 
#' @examples
#' n <- 100
#' p <- 5
#' X <- matrix(stats::rnorm(n*p),nrow=n,ncol=p)
#' #y <- stats::rbinom(n,size=1,prob=0.5)
#' y <- stats::rnorm(n)
#' glm <- glm(y~X,family="gaussian")
#' res <- residuals.glm(glm)
#' y_hat <- predict(glm,type="response")
#' all.equal(res,y-y_hat)
#' 
.residuals <- function(y,y_hat,family){
  if(length(y_hat)==1){y_hat <- matrix(y_hat,nrow=length(y),ncol=1)}
  if(family=="gaussian"){
    res <- apply(y_hat,2,function(x) (x-y)^2)
  } else if(family=="binomial"){
    if(any(!y %in% c(0,1))){stop("y outside range")}
    if(any(y_hat<0)){stop("y_hat below range")}
    if(any(y_hat>1)){stop("y_hat above range")}
    limit <- 1e-5
    y_hat[y_hat < limit] <- limit
    y_hat[y_hat > 1 - limit] <- 1 - limit
    res <- apply(y_hat,2,function(x) -2*(y*log(x)+(1-y)*log(1-x)))
  } else {
    stop("Family not implemented.")
  }
  return(res)
}

#' @title 
#' Sign discovery
#' 
#' @description
#' Assigns signs to prior weights to obtain prior coefficients
#' 
#' @inheritParams transreg
#' 
.signdisc <- function(y,X,prior,family,foldid=NULL,nfolds=10,track=FALSE){
  cond <- apply(prior,2,function(x) any(x>0) & all(x>=0))
  if(any(cond)){
    if(track){message(paste("Sign discovery procedure for source(s):",paste(which(cond),collapse=", ")))}
    #--- regression of target on features (trial) ---
    #init <- glmnet::cv.glmnet(x=X,y=y,family=family,alpha=0,foldid=foldid,nfolds=nfolds,penalty.factors=1/abs(prior))
    #sign <- sign(coef(init,s="lambda.min")[-1])
    #--- correlation between target and features (original) ---
    sign <- as.numeric(sign(stats::cor(y,X,method="spearman")))
    sign[is.na(sign)] <- 0
    #--- regression of target on prior and features (trial) ---
    #X_hat <- X * matrix(prior,nrow=nrow(X),ncol=ncol(X),byrow=TRUE)
    #object <- glmnet::glmnet(y=y,x=X_hat,lower.limits=-1,upper.limits=1,alpha=0,lambda=c(9e99,0))
    #sign <- sign(coef(object,s=0)[-1])
    # replacement
    prior[,cond] <- prior[,cond]*sign
  }
  return(prior)
}

# @example (if function 'signdisc' gets exported)
# n <- 100; p <- 500
# X <- matrix(stats::rnorm(n*p),nrow=n,ncol=p)
# beta <- stats::rnorm(p)*stats::rbinom(n=p,size=1,prob=0.2)
# y <- X %*% beta
# prior <- matrix(abs(beta),ncol=1)
# #temp <- .signdisc(y,X,prior,family="gaussian")
# #table(sign(beta),sign(temp))

#' @describeIn extract called by `coef.transreg`, `predict.transreg` and `weights.transreg`
.which.stack <- function(object,stack){
  names <- c("sta"[!is.null(object$meta.sta)],"sim"[!is.null(object$meta.sim)])
  if(is.null(stack) & length(names)==1){
    return(names)
  }
  if(!is.null(stack) & length(stack)==1 & any(names==stack)){
    return(stack)
  }
  names <- paste(paste0("stack='",names,"'"),collapse=" and ")
  stop(paste0("Choose from ",names,"."))
}

.folds <- function(y,nfolds=NULL,nfolds.ext=NULL,nfolds.int=NULL,foldid=NULL,foldid.ext=NULL,foldid.int=NULL){
  if(is.null(nfolds.int)!=is.null(nfolds.ext)|is.null(nfolds.ext)==is.null(nfolds)){
    stop("Provide either nfolds.ext and nfolds.int or nfolds.")
  }
  if(!is.null(foldid) && !all.equal(seq_len(nfolds),sort(unique(foldid)))){
    stop("nfolds and foldid are not compatible")
  }
  if(!is.null(foldid.ext) && nfolds.ext!=max(foldid.ext)){
    stop("nfolds.ext and foldid.ext are not compatible")
  }
  if(!is.null(foldid.int) && nfolds.int!=max(foldid.int)){
    stop("nfolds.int and foldid.int are not compatible")
  }
  nfolds <- list(all=nfolds,ext=nfolds.ext,int=nfolds.int)
  foldid <- list(all=foldid,ext=foldid.ext,int=foldid.int)
  for(i in c("all","ext","int")){
    if(is.null(nfolds[[i]])|!is.null(foldid[[i]])){next}
    if(all(y %in% c(0,1))){
      foldid[[i]] <- rep(x=NA,times=length(y))
      for(j in c(0,1)){
        if(i!="int"){
          foldid[[i]][y==j] <- rep(x=seq_len(nfolds[[i]]),length.out=sum(y==j))
        } else {
          if(nfolds.ext==1){
            foldid[[i]][foldid$ext==1] <- NA
            foldid[[i]][foldid$ext==0 & y==j] <- rep(x=seq_len(nfolds[[i]]),length.out=sum(foldid$ext==0 & y==j))
          } else {
            quotient <- floor(sum(y==j)/nfolds[[i]])
            remainder <- sum(y==j)%%nfolds[[i]]
            foldid[[i]][y==j] <- rep(seq_len(nfolds[[i]]),times=rep(c(quotient+1,quotient),times=c(remainder,nfolds[[i]]-remainder)))
          }
        }
      }
    } else {
      if(i!="int"){
        foldid[[i]] <- rep(x=seq_len(nfolds[[i]]),length.out=length(y))
      } else {
        if(nfolds.ext==1){
          foldid[[i]][foldid$ext==1] <- NA
          foldid[[i]][foldid$ext==0] <- rep(x=seq_len(nfolds[[i]]),length.out=sum(foldid$ext==0))
        } else {
          quotient <- floor(length(y)/nfolds[[i]])
          remainder <- length(y)%%nfolds[[i]]
          foldid[[i]] <- rep(seq_len(nfolds[[i]]),times=rep(c(quotient+1,quotient),times=c(remainder,nfolds[[i]]-remainder)))
        }
      }
    }
  }
  # permute
  if(all(y %in% c(0,1))){
    for(j in c(0,1)){
      order <- sample(seq_len(sum(y==j)))
      for(i in c("all","ext","int")){
        foldid[[i]][y==j] <- foldid[[i]][y==j][order]
      }
    }
  } else {
    order <- sample(seq_along(y))
    foldid <- lapply(foldid,function(x) x[order])
  }
  return(list(foldid=foldid[["all"]],foldid.ext=foldid[["ext"]],foldid.int=foldid[["int"]]))
}
