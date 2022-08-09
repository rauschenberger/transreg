
#'@export
#'@title
#'multiridge wrapper
#'
#'@description
#'see multiridge package.
#'
#'@param X list of matrices with n rows
#'@param Y vector of length n
#'@param family character "gaussian" or "binomial"
#'@param object multiridge-object
#'@param newx list of matrices (new data)
#'@param ... (not applicable)
#'
#'@examples
#'# simulation
#'n0 <- 100 # training samples
#'n1 <- 10000 # testing samples
#'n <- n0 + n1
#'p1 <- 5 # first covariate set
#'p2 <- 500 # second covariate set
#'X1 <- matrix(rnorm(n*p1),nrow=n,ncol=p1)
#'X2 <- matrix(rnorm(n*p2),nrow=n,ncol=p2)
#'beta1 <- rep(c(0,10),times=c(p1-4,4))
#'beta2 <- c(rnorm(100),rep(0,times=p2-100))
#'eta <- X2 %*% beta2 + X1 %*% beta1
#'family <- "binomial"
#'if(family=="gaussian"){
#'  y <- eta
#'} else if(family=="binomial"){
#'  y <- round(1/(1+exp(-eta)))
#'}
#'fold <- rep(c(0,1),times=c(n0,n1))
#'
#'# single penalty
#'glmnet <- glmnet::cv.glmnet(x=cbind(X1[fold==0,],X2[fold==0,]),y=y[fold==0],family=family,alpha=0)
#'y_hat0 <- predict(glmnet,newx=cbind(X1[fold==1,],X2[fold==1,]),s="lambda.min",type="response")
#'
#'# multiple penalties
#'object <- multiridge(X=list(X1[fold==0,],X2[fold==0,]),Y=y[fold==0],family=family)
#'y_hat1 <- predict(object,newx=list(X1[fold==1,],X2[fold==1,]))
#'
#'# comparison
#'if(family=="gaussian"){
#' loss0 <- mean((y[fold==1]-y_hat0)^2)
#' loss1 <- mean((y[fold==1]-y_hat1)^2)
#'} else if(family=="binomial"){
#' loss0 <- mean(y[fold==1]!=round(y_hat0))
#' loss1 <- mean(y[fold==1]!=round(y_hat1))
#'}
#'loss0
#'loss1
#'
#'# equivalence
#'beta <- coef(object)
#'eta2 <- beta[[1]] + X1[fold==1,] %*% beta[[2]] + X2[fold==1,] %*% beta[[3]]
#'if(family=="gaussian"){
#' y_hat2 <- eta2
#'} else if(family=="binomial"){
#' y_hat2 <- 1/(1 + exp(-eta2))
#'}
#'all.equal(y_hat1,y_hat2)
#'
#'@references
#'Mark A. van de Wiel [![MvdW](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0003-4780-8472),
#'Mirrelijn M. van Nee [![MvN](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-7715-1446),
#'and Armin Rauschenberger [![AR](https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png)](https://orcid.org/0000-0001-6498-4801) (2021).
#'"Fast Cross-validation for Multi-penalty High-dimensional Ridge Regression"
#'\emph{Journal of Computational and Graphical Statistics} 30(4):835-847
#'\url{https://doi.org/10.1080/10618600.2021.1904962}
#'
multiridge <- function(X,Y,family){
  XXblocks <- multiridge::createXXblocks(datablocks=X)
  invisible(utils::capture.output(init <- multiridge::fastCV2(XXblocks=XXblocks,Y=Y)))
  folds <- multiridge::CVfolds(Y=Y)
  invisible(utils::capture.output(final <- multiridge::optLambdasWrap(penaltiesinit=init$lambdas,
                                      XXblocks=XXblocks,Y=Y,folds=folds)))
  XXT <- multiridge::SigmaFromBlocks(XXblocks=XXblocks,penalties=final$optpen)
  object <- multiridge::IWLSridge(XXT=XXT,Y=Y,model=ifelse(family=="gaussian","linear",ifelse(family=="binomial","logistic",NA)))
  object$family <- family
  object$penalties <- final$optpen
  object$datablocks <- X
  class(object) <- "multiridge"
  return(object)
}

#'@export
#'@rdname multiridge
predict.multiridge <- function(object,newx,...){
  XXblocks <- multiridge::createXXblocks(datablocks=object$datablocks,datablocksnew=newx)
  Sigmanew <- multiridge::SigmaFromBlocks(XXblocks=XXblocks,penalties=object$penalties)
  eta <- multiridge::predictIWLS(IWLSfit=object,Sigmanew=Sigmanew)
  y_hat <- starnet:::.mean.function(eta,family=object$family)
  return(y_hat)
}

#'@export
#'@rdname multiridge
coef.multiridge <- function(object,...){
  Xblocks <- multiridge::createXblocks(datablocks=object$datablocks)
  beta <- multiridge::betasout(object,Xblocks=Xblocks,penalties=object$penalties)
  return(beta)
}
