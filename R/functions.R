
#' @export
#' 
#' @title
#' Penalised regression with multiple sets of prior coefficients
#' 
#' @description
#' Implements penalised regression with multiple sets of prior coefficients 
#' 
#' @param y
#' target: vector of length \eqn{n}
#' @param X
#' features: matrix with \eqn{n} rows and \eqn{p} columns
#' @param prior
#' prior coefficients: matrix with \eqn{p} rows and \eqn{k} columns
#' @param family
#' character "gaussian", "binomial", or "poisson"
#' @param alpha
#' elastic net mixing parameter (0=ridge, 1=lasso)
#' @param foldid
#' fold identifiers: vector of length \eqn{n} with entries from 1 to \code{nfolds}
#' @param nfolds
#' number of folds: integer
#' @param scale
#' "exp" for exponential scaling (3 parameters),
#' "iso" for isotonic scaling (1+p parameters)
#' @param sign
#' sign discovery procedure: logical
#' @param switch
#' choose between positive and negative weights for each source: logical
#' @param select
#' select from sources: logical
#' 
#' @seealso
#' Methods for objects of class \code{cornet}
#' include ...
#' 
#' @examples
#' n <- 100; p <- 500
#' X <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
#' beta <- rnorm(p)*rbinom(n=p,size=1,prob=0.2)
#' y <- X %*% beta
#' 
transreg <- function(y,X,prior,family="gaussian",alpha=1,foldid=NULL,nfolds=10,scale="iso",sign=FALSE,switch=TRUE,select=TRUE){
  
  # family <- "gaussian"; alpha <- 1; foldid <- NULL; nfolds <- 10; scale <- "exp"; sign <- FALSE; switch <- TRUE; select <- TRUE
  
  if(FALSE){
    if(!exists("family")){family <- "gaussian"}
    if(!exists("alpha")){alpha <- 1}
    if(!exists("foldid")){foldid <- NULL}
    if(!exists("nfolds")){nfolds <- 10}
    if(!exists("scale")){scale <- "iso"}
    if(!exists("sign")){sign <- FALSE}
    if(!exists("switch")){switch <- TRUE}
    if(!exists("select")){select <- TRUE}
  }
  
  if(is.vector(prior)){
    prior <- matrix(prior,ncol=1)
  }
  prior[is.na(prior)] <- 0 # allows for missing prior coefficients
  if(all(y %in% c(0,1)) != (family=="binomial")){stop("Check 'family' of 'y'.")}
  if(length(y)!=nrow(X)){stop("Entries in 'y' must match rows in 'X'.")}
  if(ncol(X)!=nrow(prior)){stop("Columns in 'X' must match rows in 'prior'.")}
  if(scale=="exp" & !switch){warning("Ignoring 'switch=FALSE'.")}
  if(all(prior>=0) & !sign){warning("Consider sign discovery procedure.")}
  
  if(is.null(foldid)){
    foldid <- palasso:::.folds(y=y,nfolds=nfolds)
  }
  
  n <- nrow(X); p <- ncol(X)
  k <- ncol(prior)
  nfolds <- max(foldid)
  
  base <- glmnet::cv.glmnet(y=y,x=X,family=family,alpha=alpha,nlambda=100,keep=TRUE,foldid=foldid)
  
  #full <- glmnet::glmnet(y=y,x=X,family=family,alpha=alpha) # trial
  ##all(coef(base,s=base$lambda.min)==coef(full,s=base$lambda.min))
  #Y_hat <- int <- matrix(NA,nrow=n,ncol=length(full$lambda)) # trial
  
  # external re-scaling
  prior.ext <- prior
  if(sign){
    prior.ext <- sign.disc(y=y,X=X,prior=prior.ext,family=family)
  }
  base$z <- prior.ext
  if(scale=="exp"){
    prior.ext <- exp.multiple(y=y,X=X,prior=prior.ext,family=family,select=select)
  } else if(scale=="iso"){
    prior.ext <- iso.multiple(y=y,X=X,prior=prior.ext,family=family,switch=switch,select=select)
  } else {
    stop("Invalid scale.",call.=FALSE)
  }
  
  # start trial nested cv
  # fit glmnet on meta-features and features, using all data
  #temp <- X %*% prior.ext$beta
  #full <- glmnet::glmnet(y=y,x=cbind(temp,X),alpha=alpha,family=family,lower.limits=rep(c(0,-Inf),times=c(k,p)),
  #                  penalty.factor=rep(c(0,1),times=c(k,p)))
  #pred <- matrix(NA,nrow=n,ncol=length(full$lambda))
  # end trial nested cv
  
  y_hat <- matrix(NA,nrow=n,ncol=k+2) # was ncol=k+1
  for(i in seq_len(nfolds)){
    y0 <- y[foldid!=i]
    X0 <- X[foldid!=i,]
    X1 <- X[foldid==i,]
    
    # internal re-scaling
    prior.int <- prior
    if(sign){
      prior.int <- sign.disc(y=y0,X=X0,prior=prior.int,family=family)
      # trial:
      #prior.int[prior.ext<0] <- -prior.int[prior.ext<0] # temporary
    }
    if(scale=="exp"){
      prior.int <- exp.multiple(y=y0,X=X0,prior=prior.int,family=family,select=select)
    } else if(scale=="iso"){
      prior.int <- iso.multiple(y=y0,X=X0,prior=prior.int,family=family,switch=switch,select=select)
    } else {
      stop("Invalid scale.",call.=FALSE)
    }
    
    # predictions
    for(j in seq_len(k)){
      y_hat[foldid==i,j] <- X1 %*% prior.int$beta[,j] # original (harmonise with predict.transreg)
      #y_hat[foldid==i,j] <- prior.int$alpha[j] + X1 %*% prior.int$beta[,j] # trial 2022-01-04 (see below)
    }
    y_hat[foldid==i,k+1] <- base$fit.preval[base$foldid==i,base$lambda==base$lambda.min]
    y_hat[foldid==i,k+2] <- base$fit.preval[base$foldid==i,base$lambda==base$lambda.1se]
    
    #part <- glmnet::glmnet(y=y0,x=X0,family=family,alpha=alpha) # trial
    #Y_hat[foldid==i,] <- stats::predict(part,s=full$lambda,newx=X1) # trial
    #int[foldid==i,] <- rep(coef(part,s=full$lambda)["(Intercept)",],each=sum(foldid==i)) # trial
    
    # start trial nested cv
    # fit glmnet on meta-features and features, using included folds
    #temp <- X0 %*% prior.int$beta
    #sub <- glmnet::glmnet(y=y0,x=cbind(temp,X0),alpha=alpha,family=family,lower.limits=rep(c(0,-Inf),times=c(k,p)),
    #                       penalty.factor=rep(c(0,1),times=c(k,p)))
    #temp <- X1 %*% prior.int$beta
    #pred[foldid==i,] <- predict(sub,newx=cbind(temp,X1),s=full$lambda)
    # end trial nested cv
  }
  
  #cvm <- palasso:::.loss(y=y,fit=joinet:::.mean.function(Y_hat,family=family),family=family,type.measure="deviance",foldid=foldid)[[1]] # trial
  #id <- which.min(cvm) # trial
  #temp <- Y_hat[,id] # trial
  #all(y_hat[,k+1]==temp) # trial
  #y_hat[,k+1] <- Y_hat[,id]-int[,id] # trial (harmonise with predict.transreg)
  
  base$prior <- prior.ext
  meta <- glmnet::cv.glmnet(y=y,x=y_hat,foldid=foldid,alpha=1,family=family,lower.limits=0)
  
  # start alternative meta-features
  trial <- glmnet::cv.glmnet(y=y,x=cbind(y_hat[,1:k],X),alpha=alpha,family=family,
                             lower.limits=rep(c(0,-Inf),times=c(k,p)),
                             penalty.factor=rep(c(0,1),times=c(k,p)),foldid=foldid)
  
  #trial <- glmnet::cv.glmnet(y=y,x=cbind(y_hat[,1:k],X),alpha=alpha,family=family,lower.limits=rep(c(0,-Inf),times=c(k,p)),
  #                           penalty.factor=rep(c(0,1),times=c(k,p)))
  # tune weight here (penalty factor = inverse weight):
  # weight between 0 and 1 for meta-features, weight 1 for original features
  #pf <- 0 #seq(from=0,to=1,length.out=6)
  #trials <- list()
  #for(i in seq_along(pf)){
  #  trials[[i]] <- glmnet::cv.glmnet(y=y,x=cbind(y_hat[,1:k],X),alpha=alpha,family=family,lower.limits=rep(c(0,-Inf),times=c(k,p)),
  #                                   penalty.factor=rep(c(pf[i],1),times=c(k,p)),foldid=foldid)
  #}
  ## try: weight <- seq(from=0.5,to=1,length.out=6)
  ## try: penalty.factor=c(1/(1-weight[i]),1/weight[i])
  #cvm <- sapply(trials,function(x) min(x$cvm))
  #cat("min cvm:",cvm)
  #tryCatch(expr=plot(x=pf,y=cvm),error=function(x) NULL)
  #test <- trials[[which(pf==0)]]
  # end alternative meta-features
  
  # start alternative tune pf
  # test <- trials[[which.min(cvm)]]
  # end alternative tune pf
  
  # start alternative no cv
  #temp <- X %*% prior.ext$beta
  #test <- glmnet::cv.glmnet(y=y,x=cbind(temp,X),alpha=alpha,family=family,lower.limits=rep(c(0,-Inf),times=c(k,p)),
  #                                 penalty.factor=rep(c(0,1),times=c(k,p)),foldid=foldid)
  # end alternative no cv
  
  # start trial nested cv
  #loss <- apply(pred,2,function(x) starnet:::.loss(y=y,x=x,family=family,type.measure="deviance"))
  #lambda.min <- full$lambda[which.min(loss)]
  #test <- list(glmnet.fit=full,lambda.min=lambda.min)
  #class(test) <- "cv.glmnet"
  # temp <- X %*% prior.ext$beta
  # predict(test,newx=cbind(temp,X),s="lambda.min")
  # end trial nested cv
  
  # start alternative triple
  #test <- glmnet::cv.glmnet(y=y,x=cbind(y_hat,X),alpha=alpha,family=family,lower.limits=rep(c(0,-Inf),times=c(k+2,p)),
  #                                 penalty.factor=rep(c(0,1),times=c(k+2,p)),foldid=foldid)
  # end alternative triple
  
  object <- list(base=base,meta=meta,scale=scale,trial=trial) # test=test
  class(object) <- "transreg"
  return(object)
}

#' @export
#' 
#' @title
#' Predict
#' 
#' @description
#' Predicts outcome
#' 
#' @param object
#' object of class 'transreg'
#' @param newx
#' features: matrix with m rows (samples) and p columns (variables)
#' @param ...
#' not applicable
#' 
#' @examples
#' NA
#' 
predict.transreg <- function(object,newx,...){
  one <- newx %*% object$base$prior$beta # original (harmonise with transreg)
  #one <- object$base$prior$alpha + newx %*% object$base$prior$beta # trial 2022-01-04 (see above)
  two <- stats::predict(object$base,s=c(object$base$lambda.min,object$base$lambda.1se),newx=newx) # was s="lambda.min"
  #two <- newx %*% coef(object$base,s="lambda.min")[-1] # trial (harmonise with transreg)
  y_hat <- stats::predict(object$meta,s="lambda.min",newx=cbind(one,two),type="response")
  return(y_hat)
}

predict.trial <- function(object,newx,...){
  one <- newx %*% object$base$prior$beta
  y_hat <- stats::predict(object$trial,s="lambda.min",newx=cbind(one,newx),type="response")
  return(y_hat)
}

#predict.test <- function(object,newx,...){
#  one <- newx %*% object$base$prior$beta
#  y_hat <- stats::predict(object$test,s="lambda.min",newx=cbind(one,newx),type="response")
#  return(y_hat)
#}

#predict.test <- function(object,newx,...){
#  one <- newx %*% object$base$prior$beta
#  two <- stats::predict(object$base,s=c(object$base$lambda.min,object$base$lambda.1se),newx=newx) # was s="lambda.min"
#  y_hat <- stats::predict(object$test,s="lambda.min",newx=cbind(one,two,newx),type="response")
#  return(y_hat)
#}

#' @export
#'
#' @title Calculate residuals
#' 
#' @description
#' Calculates residuals from observed outcome
#' and predicted values (Gaussian family)
#' or predicted probablities (binomial family).
#' 
#' @inheritParams transreg
#' @param y_hat
#' predicted values: vector of length n, or matrix with n rows and k columns
#' @param family
#' "gaussian" or "binomial"
#' 
#' @examples
#' NA
#' 
residuals <- function(y,y_hat,family){
  if(length(y_hat)==1){y_hat <- matrix(y_hat,nrow=length(y),ncol=1)}
  if(family=="gaussian"){
    res <- apply(y_hat,2,function(x) (x-y)^2)
  } else if(family=="binomial"){
    if(any(y_hat<0)){stop("Below range.")}
    if(any(y_hat>1)){stop("Above range.")}
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
#' @examples
#' NA
#' 
sign.disc <- function(y,X,prior,family,foldid=NULL,nfolds=10){
  cond <- apply(prior,2,function(x) any(x>0) & all(x>=0))
  if(any(cond)){
    message(paste("Sign discovery procedure for source(s):",paste(which(cond),collapse=", ")))
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

#' @title 
#' 
#' Exponential scaling
#' 
#' @description
#' Performs exponential scaling
#' 
#' @inheritParams transreg
#' @param plot logical
#' 
#' @examples
#' NA
#' 
exp.multiple <- function(y,X,prior,family,select,plot=TRUE){
  
  #message("Exponential scaling ...")
  
  n <- nrow(X); p <- ncol(X); k <- ncol(prior)
  
  alpha <- gamma <- pvalue <- sign <- tau <- rep(NA,times=k)
  beta <- matrix(NA,nrow=p,ncol=k)
  
  # This should be done with training data, within cross-validation loop, for each set of prior information. Also consider more complex transformations.
  
  exp <- c(0,exp(seq(from=log(0.01),to=log(50),length.out=100)))
  
  for(j in seq_len(k)){
    coefs <- matrix(NA,nrow=length(exp),ncol=2)
    pred <- matrix(NA,nrow=nrow(X),ncol=length(exp))
    for(i in seq_along(exp)){
      temp <- sign(prior[,j])*abs(prior[,j])^exp[i]
      eta <- X %*% temp
      # Use simple linear/logistic/Poisson regression of y on eta, and extract fitted values. This should solve the scaling issue. => But be aware of negative coefficients!
      glm <- stats::glm(y~eta,family=family)
      #coefs[i,] <- stats::coef(glm) # original
      temp <- stats::coef(glm)
      coefs[i,1] <- temp["(Intercept)"] # trial
      coefs[i,2] <- ifelse(is.na(temp["eta"]),0,temp["eta"]) # trial
      pred[,i] <- stats::fitted(glm)
    }
    cvm <- palasso:::.loss(y=y,fit=pred,family=family,type.measure="deviance")[[1]]
    id.min <- which.min(cvm)
    
    alpha[j] <- coefs[id.min,1]
    gamma[j] <- coefs[id.min,2]
    beta[,j] <- coefs[id.min,2]*sign(prior[,j])*abs(prior[,j])^exp[id.min]
    sign[j] <- sign(coefs[id.min,2])
    tau[j] <- exp[id.min]
    
    # This was going completely wrong:
    #a <- pred[,which.min(cvm)]
    #b <- X %*% prior[,j]
    #plot(x=a,y=b)
    # These things should be perfectly positively correlated!
    
    # Include here test of whether residuals are significantly lower than those of the empty model.
    if(select){
      res0 <- residuals(y=y,y_hat=mean(y),family=family)
      res1 <- residuals(y=y,y_hat=pred[,which.min(cvm),drop=FALSE],family=family)
      pvalue[j] <- stats::wilcox.test(x=res1,y=res0,paired=TRUE,alternative="less")$p.value
    }
    
    if(plot){
      tryCatch(graphics::plot(x=exp,y=cvm,main=j),error=function(x) NULL)
      tryCatch(graphics::abline(v=exp[id.min]),error=function(x) NULL)
    }
  }
  
  if(select){
    remove <- pvalue > 0.05/k
    # Use same cut-off for exponential and isotonic scaling!
    beta[,remove] <- 0
    message(ifelse(remove,".",ifelse(sign==1,"+","-")))
  }
  
  return(list(alpha=alpha,beta=beta,gamma=gamma,tau=tau))
}

#' @title 
#' Isotonic scaling
#' 
#' @description 
#' Performs isotonic scaling. This function is for comparison only.
#'
#' @inheritParams transreg
#' 
#' @seealso
#' iso.fast.single
#' 
#' @examples 
#' NA
#' 
iso.slow.single <- function(y,X,prior,family){
  
  message("Isotonic scaling ...")
  
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
      BETA[,i] <- 0 #was prior[,i]
    } else {
      ALPHA[i] <- result$getValue(alpha)
      beta.r <- result$getValue(beta)
      BETA[,i] <- beta.r[order(order)]
    }
    #post <- stats::coef(stats::glm(y~X,family="binomial"))
    #cbind(post,c(alpha.r,posterior[,i]))
  }
  
  BETA[-1e-06 < BETA & BETA < +1e-06] <- 0
  return(list(alpha=ALPHA,beta=BETA))
}

# with sign constraints (glmnet)
#' @title 
#' Isotonic scaling
#' 
#' @description 
#' Performs isotonic scaling
#' 
#' @inheritParams transreg
#' 
#' @examples
#' NA
#' 
iso.fast.single <- function(y,X,prior,family){
  n <- length(y)
  p <- nrow(prior)
  k <- ncol(prior)
  ALPHA <- rep(NA,times=k)
  BETA <- posterior <- matrix(data=NA,nrow=p,ncol=k)
  for(i in seq_len(k)){
    single.prior <- prior[,i]
    order <- order(single.prior)
    Xo <- X[,order,drop=FALSE]
    #cond <- single.prior[order]>=0
    sign <- sign(single.prior[order])
    A <- t(apply(Xo[,sign<0,drop=FALSE],1,function(x) cumsum(x))) # was !cond
    B <- t(apply(Xo[,sign>=0,drop=FALSE],1,function(x)  rev(cumsum(rev(x))))) # was cond
    if(any(sign<0)&any(sign>=0)){ # new condition
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
    #Xcs <- cbind(A,B) # original
    times <- c(sum(sign==-1),sum(sign==0),sum(sign==+1)) # was c(sum(!cond),sum(cond)) 
    # Bug fix on 2022-03-10: replaced family="gaussian" by family=family
    inc <- glmnet::glmnet(y=y,x=Xcs,family=family,alpha=0,lambda=0,lower.limits=rep(c(-Inf,0,0),times=times),upper.limits=rep(c(0,0,Inf),times=times))
    coef <- stats::coef(inc)
    ALPHA[i] <- coef[1]
    beta <- coef[-1]
    new <- c(rev(cumsum(rev(beta[sign<0]))),cumsum(beta[sign>=0])) # was !cond and cond
    BETA[,i] <- new[order(order)]
  }
  return(list(alpha=ALPHA,beta=BETA))
}

# Speed up this function by only examining the "negative prior" if the "positive prior" is insignificant? Or do some pre-screening based on correlation, and then decide which one to run first and which one to run second (i.e. only if the first one fits poorly).

#' @title 
#' Isotonic scaling
#' 
#' @description
#' Performs isotonic scaling
#' 
#' @inheritParams transreg
#' 
#' @examples 
#' NA
#' 
iso.multiple <- function(y,X,prior,family,switch=TRUE,select=TRUE){
  
  k <- ncol(prior)
  
  #ALPHA <- rep(NA,times=ncol(prior))
  #BETA <- matrix(NA,nrow=nrow(prior),ncol=ncol(prior))
  
  prior0 <- iso.fast.single(y=y,X=X,prior=+prior,family=family) # was iso.slow.single
  alpha0 <- prior0$alpha; beta0 <- prior0$beta
  # beta0[beta0==min(beta0)] <- sort(beta0,decreasing=FALSE)[2] # hack (trial)
  fit0 <- joinet:::.mean.function(alpha0 + X %*% beta0,family=family)
  
  res <- residuals(y=y,y_hat=mean(y),family=family)
  res0 <- residuals(y=y,y_hat=fit0,family=family)
  
  #if(family=="gaussian"){
  #  res <- (y-mean(y))^2
  #  res0 <- apply(fit0,2,function(x) (x-y)^2)
  #} else if(family=="binomial"){
  #  res <- -2*(y*log(mean(y))+(1-y)*log(1-mean(y)))
  #  #if(any(fit0>1)|any(fit0<0)|any(fit1>1)|any(fit1<0)){stop("invalid.")}
  #  limit <- 1e-5
  #  fit0[fit0 < limit] <- limit
  #  fit0[fit0 > 1 - limit] <- 1 - limit
  #  res0 <- apply(fit0,2,function(x) -2*(y*log(x)+(1-y)*log(1-x)))
  #} else {
  #  stop("Family not implemented.")
  #}
  #loss0 <- apply(fit0,2,function(x) palasso:::.loss(y=y,fit=x,family=family,type.measure="deviance")[[1]])
  # when writing a function for residuals: check loss0==colMeans(res0)
  
  pval0 <- apply(res0,2,function(x) stats::wilcox.test(x=x,y=res,paired=TRUE,alternative="less")$p.value)
  
  if(switch){
    prior1 <- iso.fast.single(y=y,X=X,prior=-prior,family=family) # was iso.slow.single
    alpha1 <- prior1$alpha; beta1 <- prior1$beta
    fit1 <- joinet:::.mean.function(alpha1 + X %*% beta1,family=family)
    
    res1 <- residuals(y=y,y_hat=fit1,family=family)
    
    #if(family=="gaussian"){
    #  res1 <- apply(fit1,2,function(x) (x-y)^2)
    #} else if(family=="binomial"){
    #  limit <- 1e-5
    #  fit1[fit1 < limit] <- limit
    #  fit1[fit1 > 1 - limit] <- 1 - limit
    #  res1 <- apply(fit1,2,function(x) -2*(y*log(x)+(1-y)*log(1-x)))
    #} else {
    #  stop("Family not implemented.")
    #}
    
    pval1 <- apply(res1,2,function(x) stats::wilcox.test(x=x,y=res,paired=TRUE,alternative="less")$p.value)
    
    cat(paste(signif(pval0,digits=2),sep=" "),"\n")
    cat(paste(signif(pval1,digits=2),sep=" "),"\n")
    
    cond <- pval0 <= pval1
    #ALPHA[cond] <- alpha0[cond]
    #ALPHA[!cond] <- alpha0[!cond]
    #prior[,cond] <- beta0[,cond] # original
    #prior[,!cond] <- beta1[,!cond] # original
    alpha0[!cond] <- alpha1[!cond] # trial 2022-01-13
    beta0[,!cond] <- beta1[,!cond] # correction
    pvalue <- pmin(pval0,pval1)
    
  } else {
    cond <- rep(TRUE,times=k)
    pvalue <- pval0
    res1 <- res0 # temporary
  }
  
  if(select){
    #remove <- pvalue>0.05/((1+switch)*k) # original
    remove <- pvalue>0.05 # trial
    
    # Alternatively, include all sources that lead to any improvement
    # (no matter how small).
    # Or remove multiple testing correction.
    #warning("Temporary line:")
    #remove <- pmin(colMeans(res0),colMeans(res1)) >= mean(res) # trial
    #prior[,remove] <- 0 # changed from 1 (mistake!?) to 0 # original
    alpha0[remove] <- 0 # correction
    beta0[,remove] <- 0 # correction
    message(ifelse(remove,".",ifelse(cond,"+","-")))
  }
  
  return(list(alpha=alpha0,beta=beta0)) # was alpha=NULL # trial 2022-01-04
}

#' @export
#' 
#' @title 
#' Cross-validation
#' 
#' @description
#' Performs external \eqn{k}-fold cross-validation
#' to estimate the predictive performance of different methods
#' 
#' @inheritParams transreg
#' @param target
#' list with slot x (feature matrix with n rows and p columns) and slot y (target vector of length n)
#' @param source
#' list of k lists, each with slot x (feature matrix with m_i rows and p columns) and slot y (target vector of length m_i)
#' @param partitions monotone: for GRridge
#' @param alpha.prior number between 0 (lasso) and 1 (ridge), character "p-value", or NULL (alpha.prior=alpha, but if alpha=1 then alpha.prior=0.95)
#' @param z prior weights
#' @param foldid.ext external fold identifiers
#' @param nfolds.ext number of external folds
#' @param foldid.int internal fold identifiers
#' @param nfolds.int number of internal folds
#' @param type.measure character
#' @param alpha.prior alpha for source regression
#' @param monotone logical
#' 
#' @examples
#' NA
#' 
cv.transfer <- function(target,source=NULL,prior=NULL,z=NULL,family,alpha,scale="iso",sign=FALSE,select=TRUE,switch=TRUE,foldid.ext=NULL,nfolds.ext=10,foldid.int=NULL,nfolds.int=10,type.measure="deviance",alpha.prior=NULL,partitions=NULL,monotone=NULL){
  
  if(FALSE){
    if(!exists("source")){source <- NULL}
    if(!exists("prior")){prior <- NULL}
    if(!exists("z")){z <- NULL}
    if(!exists("foldid.ext")){foldid.ext <- NULL}
    if(!exists("nfolds.ext")){nfolds.ext <- 10}
    if(!exists("foldid.int")){foldid.int <- NULL}
    if(!exists("nfolds.int")){nfolds.int <- 10}
    if(!exists("type.measure")){type.measure <- "deviance"}
    if(!exists("alpha.prior")){alpha.prior <- NULL}
    if(!exists("partitions")){partitions <- NULL}
    if(!exists("monotone")){monotone <- NULL}
  }
  
  if(is.null(source)==is.null(prior)){
    stop("Provide either \"source\" or \"prior\".",call.=FALSE)
  }
  
  if(!is.list(target)||is.null(names(target))){stop("Expect argument target as list with slots x and y.")}
  #names(target) <- tolower(names(target))
  if(!any(names(target) %in% c("y","x"))){stop("Expect argument target as list with slots x and y.")}
  
  if(!is.null(source)){
  if(!is.list(source)){stop("Expect argument source as list of lists.")}
  for(i in seq_along(source)){
    if(!is.list(source[[i]])||is.null(names(source[[i]]))){stop("Expect argument source as list of lists with slots x and y.")}
    #names(source[[i]]) <- tolower(names(source[[i]]))
    if(!any(names(source[[i]]) %in% c("y","x"))){stop("Expect argument source as list of lists with slots x and y.")}
  }
  }
  
  if(is.vector(prior)){
    prior <- matrix(prior,ncol=1)
  }
  
  if(!exists(".Random.seed")){
    .Random.seed <- NULL
  }
  seed <- .Random.seed
  
  #  family <- "binomial"; alpha <- 1; foldid.int <- foldid.ext <- NULL; nfolds.ext <- 2; nfolds.int <- 10; type.measure <- "deviance"; scale <- "iso"; sign <- FALSE; select <- FALSE; alpha.prior <- NULL
  
  #n <- c(nrow(target$x),vapply(source,function(x) nrow(x$x),numeric(1)))
  #N <- matrix(n,nrow=nrow(sds),ncol=ncol(sds))
  #mus <- cbind(colMeans(target$x),vapply(source,function(x) colMeans(x$x),numeric(ncol(target$x))))
  #sds <- cbind(apply(target$x,2,sd),vapply(source,function(x) apply(x$x,2,sd),numeric(ncol(target$x))))
  #sd <- sqrt(rowSums(sds^2 * N)/n)
  #mu <- rowSums(mus)/n
  
  if(is.data.frame(target$x)){
    target$x <- as.matrix(target$x)
  }
  #if(family=="gaussian"){target$y <- as.numeric(scale(target$y))} # temporary
  
  if(!is.null(source)){
    standardise <- "old" # Changed from "old" to "new" on 2022-01-06.
    
    if(standardise=="old"){
      target$x <- scale(target$x) # original
      target$x[is.na(target$x)] <- 0 # original
    }
    
    for(i in seq_along(source)){
      if(is.data.frame(source[[i]]$x)){
        source[[i]]$x <- as.matrix(source[[i]]$x)
      }
      #if(family=="gaussian"){source[[i]]$y <- as.numeric(scale(source[[i]]$y))} # temporary
      if(standardise=="old"){
        source[[i]]$x <- scale(source[[i]]$x) # original
        source[[i]]$x[is.na(source[[i]]$x)] <- 0 # original
      }
    }
    
    # CONTINUE HERE: standardise target x and all source x's with overall mean and overall variance (across target x and all source x's).
    if(standardise=="new"){
      temp <- rbind(target$x,do.call(what="rbind",args=lapply(source,function(x) x$x)))
      mu_temp <- colMeans(temp)
      sd_temp <- apply(temp,2,sd)
      n <- c(length(target$y),sapply(source,function(x) length(x$y)))
      mus <- rbind(colMeans(target$x),t(sapply(source,function(x) colMeans(x$x))))
      #sds <- rbind(apply(target$x,2,sd),t(sapply(source,function(x) apply(x$x,2,sd))))
      mu <- colSums(mus*n)/sum(n)
      sd <- sqrt(rowSums(cbind(rowSums((t(target$x)-mu)^2),sapply(source,function(x) rowSums((t(x$x)-mu)^2))))/(sum(n)-1))
      if(any(abs(mu-mu_temp)>1e-06)){stop("Invalid.")}
      if(any(abs(sd-sd_temp)>1e-06)){stop("Invalid.")}
      
      target$x <- t((t(target$x)-mu)/sd)
      target$x[,sd==0] <- 0
      for(i in seq_along(source)){
        source[[i]]$x <- t((t(source[[i]]$x)-mu)/sd)
        source[[i]]$x[,sd==0] <- 0
      }
    }
    
    if(is.null(alpha.prior)){alpha.prior <- ifelse(alpha==1,0.95,alpha)}
    
    # prior
    k <- length(source)
    p <- ncol(target$x)
    prior <- matrix(NA,nrow=p,ncol=k)
    for(i in seq_len(k)){
      foldid.source <- palasso:::.folds(source[[i]]$y,nfolds=10)
      #temp <- scale(source[[i]]$x)
      #temp[is.na(temp)] <- 0
      if(is.character(alpha.prior) & alpha.prior=="p-value"){
        if(family=="binomial"){
          p.value <- apply(source[[i]]$x,2,function(x) stats::wilcox.test(x~source[[i]]$y)$p.value)
          sign <- apply(source[[i]]$x,2,function(x) mean(x[source[[i]]$y==1])-mean(x[source[[i]]$y==0]))
          prior[,i] <- sign*(-log10(p.value))
          #test <- cor(source[[i]]$y,source[[i]]$x)
        } else {
          stop("Implement correlation test!",call.=FALSE)
        }
      } else {
        object <- glmnet::cv.glmnet(y=source[[i]]$y,x=source[[i]]$x,family=family,alpha=alpha.prior,foldid=foldid.source) # was alpha=ifelse(alpha==1,0.95,alpha) # added scale()
        # Why not always ridge regression for prior?
        prior[,i] <- as.numeric(stats::coef(object,s="lambda.min"))[-1]
      }
    }
    
  } else {
    k <- ncol(prior)
  } # end if(!is.null(source))
  
  if(!is.null(partitions)){trial <- TRUE}else{trial <- FALSE}
  if(!is.null(z)){trial2 <- TRUE}else{trial2 <- FALSE}
  
  #if(trial){
  #  warning("Temporary: absolute values (for testing purposes)")
  #  prior <- abs(prior) # remove this!!!
  #}
  
  #--- fold identifiers ---
  if(is.null(foldid.ext)){
    foldid.ext <- palasso:::.folds(y=target$y,nfolds=nfolds.ext)
  } else {
    nfolds.ext <- max(foldid.ext)
  }
  
  names <- c("mean","glmnet","glmtrans"[!is.null(source)],"transreg","transreg.trial","GRridge"[trial],"NoGroups"[trial],"fwelnet"[trial2],"xtune"[trial2],"CoRF"[trial2],"ecpc"[trial2]) # "transreg.test"
  pred <- matrix(data=NA,nrow=length(target$y),ncol=length(names),dimnames=list(NULL,names))
  time <- rep(0,time=length(names)); names(time) <- names
  
  set.seed(seed) # trial
  for(i in seq_len(nfolds.ext)){
    y0 <- target$y[foldid.ext!=i]
    X0 <- target$x[foldid.ext!=i,]
    X1 <- target$x[foldid.ext==i,]
    if(is.null(foldid.int)){
      foldid <- palasso:::.folds(y=y0,nfolds=nfolds.int)
    } else {
      foldid <- foldid.int[foldid.ext!=i]
    }
    
    # mean
    set.seed(seed)
    start <- Sys.time()
    pred[foldid.ext==i,"mean"] <- rep(mean(y0),times=sum(foldid.ext==i))
    end <- Sys.time()
    time["mean"] <- time["mean"]+difftime(time1=end,time2=end,units="secs")
    
    # glmnet
    set.seed(seed) # trial
    start <- Sys.time()
    object <- glmnet::cv.glmnet(y=y0,x=X0,family=family,foldid=foldid,alpha=alpha)
    pred[foldid.ext==i,"glmnet"] <- as.numeric(stats::predict(object,newx=X1,s="lambda.min",type="response"))
    end <- Sys.time()
    time["glmnet"] <-  time["glmnet"]+difftime(time1=end,time2=start,units="secs")
    
    # glmtrans
    if(!is.null(source)){
      set.seed(seed) # trial
      start <- Sys.time()
      object <- tryCatch(glmtrans::glmtrans(target=list(x=X0,y=y0),source=source,alpha=alpha,family=family,nfolds=nfolds.int),error=function(x) NULL)
      if(!is.null(object)){
        pred[foldid.ext==i,"glmtrans"] <- stats::predict(object,newx=X1,type="response")
      }
      end <- Sys.time()
      time["glmtrans"] <-  time["glmtrans"]+difftime(time1=end,time2=start,units="secs")
    }
    
    # transreg
    set.seed(seed) # trial
    start <- Sys.time()
    object <- transreg(y=y0,X=X0,prior=prior,family=family,foldid=foldid,alpha=alpha,scale=scale,sign=sign,switch=switch,select=select)
    pred[foldid.ext==i,"transreg"] <- stats::predict(object,newx=X1)
    end <- Sys.time()
    time["transreg"] <- time["transreg"]+difftime(time1=end,time2=start,units="secs")
    
    # transreg trial
    pred[foldid.ext==i,"transreg.trial"] <- predict.trial(object,newx=X1)
    #pred[foldid.ext==i,"transreg.test"] <- predict.test(object,newx=X1)
    
    # GRridge
    if(trial){
      set.seed(seed) # trial
      start <- Sys.time()
      #index <- list()
      #for(j in seq_len(ncol(prior))){
      #  index[[j]] <- GRridge::CreatePartition(vec=abs(prior[,j]),ngroup=10,mingr=pmin(25,nrow(prior)/20))
      #}
      # Why do the three grouping have groups of exactly the same size? => Mistake!
      #monotone <- rep(TRUE,times=length(index))
      object <- tryCatch(GRridge::grridge(highdimdata=t(X0),response=y0,
                                          partitions=partitions,monotone=monotone,
                                          innfold=nfolds.int,fixedfoldsinn=TRUE),
                         error=function(x) NULL)
      if(!is.null(object)){
        temp <- GRridge::predict.grridge(object=object,datanew=as.data.frame(t(X1)))
        pred[foldid.ext==i,"NoGroups"] <- temp[,"NoGroups"]
        pred[foldid.ext==i,"GRridge"] <- temp[,"GroupRegul"]
      }
      
      end <- Sys.time()
      time["GRridge"] <- time["GRridge"]+difftime(time1=end,time2=start,units="secs")
    }
    
    # co-data methods
    if(trial2){
      
      # fwelnet
      set.seed(seed)
      start <- Sys.time()
      #if(is.vector(z)){z <- as.matrix(z,ncol=1)}
      object <- tryCatch(fwelnet::cv.fwelnet(x=X0,y=y0,z=z,family=family,foldid=foldid,alpha=alpha),error=function(x) NULL)
      if(!is.null(object)){
        pred[foldid.ext==i,"fwelnet"] <- stats::predict(object,xnew=X1,s="lambda.min",type="response")
      }
      end <- Sys.time()
      time["fwelnet"] <- time["fwelnet"]+difftime(time1=end,time2=start,units="secs")

      # ecpc
      set.seed(seed)
      start <- Sys.time()
      if(is.vector(z)){
        Z <- list(as.numeric(z))
      } else {
        Z <- as.list(as.data.frame(z))
      }
      object <- tryCatch(ecpc::ecpc(Y=y0,X=X0,Z=Z,X2=X1,fold=nfolds.int),error=function(x) NULL)
      if(!is.null(object)){
        pred[foldid.ext==i,"ecpc"] <- object$Ypred
      }
      end <- Sys.time()
      time["ecpc"] <- time["ecpc"]+difftime(time1=end,time2=start,units="secs")
      
      # # xtune
      # set.seed(seed)
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

      # # CoRF
      # set.seed(seed)
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
  
  #loss <- palasso:::.loss(y=target$y[foldid.ext!=0],fit=pred[foldid.ext!=0,],
  #                        family=family,type.measure=type.measure)[[1]]
  
  loss <- apply(pred,2,function(x) starnet:::.loss(y=target$y[foldid.ext!=0],x=x[foldid.ext!=0],
                          family=family,type.measure=type.measure))
  
  attributes(loss)$time <- time
  
  return(loss)
}

# CONTINUE HERE: cv.transfer function (so far only used for single split, error if other fold identifiers than 0 and 1?)

#' @title
#' Simulation
#' 
#' @description 
#' Simulates data
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
#' 
#' @examples
#' NA
#' 
.simulate <- function(p=1000,n.target=100,n.source=150,k=3,family="gaussian",prop=0.01,rho.beta=0.95,rho.x=0.95,w=0.5){
  
  target <- source <- list()
  
  mu <- rep(x=0,times=k+1)
  Sigma <- matrix(data=rho.beta,nrow=k+1,ncol=k+1) # original
  diag(Sigma) <- 1 # original
  #Sigma <- matrix(data=NA,nrow=k+1,ncol=k+1)
  #Sigma <- rho.beta^(abs(row(Sigma)-col(Sigma))) # trial 2022-01-10
  beta <- mvtnorm::rmvnorm(n=p,mean=mu,sigma=Sigma)
  
  if(FALSE){
    message("Temporary re-scaling of coefficients!")
    
    # source (with perturbation: exponentiated effects)
    #exp <- 0.2
    #beta[,1] <- sign(beta[,1])*abs(beta[,1])^exp
    #beta[,2] <- sign(beta[,2])*abs(beta[,2])^exp
    #beta[,3] <- sign(beta[,3])*abs(beta[,3])^exp
    
    # source (with perturbation: partially non-informative source data)
    #beta[,1] <- sample(beta[,1])
    #beta[,2] <- sample(beta[,2])
    #beta[,3] <- beta[,3]
    
    # target
    #beta[,4] <- beta[,4] # should not change
  }
  
  if(FALSE){
    warning("Remove the following line:")
    Sigma <- 0*Sigma
    diag(Sigma) <- 1
  }
  
  cond <- mvtnorm::rmvnorm(n=p,mean=mu,sigma=Sigma)>stats::qnorm(1-prop)
  beta[cond==0] <- 0
  
  #stats::cor(beta,method="spearman")
  
  # features
  #w <- 0.50 # weight between signal and noise
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
  
  #foldid.ext <- rep(c(0,1),times=c(n0,n1))
  
  return(list(source=source,target=target,beta=beta))
}
