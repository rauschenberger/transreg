
set.seed(1)
n <- 100; p <- 1000
X <- matrix(stats::rnorm(n*p),nrow=n,ncol=p)
beta <- stats::rnorm(p)
y <- X %*% beta

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#- - - slow and fast isotonic scaling- - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

prior1 <- beta + stats::rnorm(p,sd=0.1)
prior2 <- beta + stats::rnorm(p,sd=0.1)

slow <- iso.slow.single(y=y,X=X,prior=matrix(prior1,ncol=1),family="gaussian")$beta
fast <- iso.fast.single(y=y,X=X,prior=matrix(prior1,ncol=1),family="gaussian")$beta

testthat::test_that("expected signs (slow)",{
  cond1 <- all(slow[prior1==0]==0)
  cond2 <- all(slow[prior1<0]<=0)
  cond3 <- all(slow[prior1>0]>=0)
  testthat::expect_true(cond1&cond2&cond3)
})

testthat::test_that("expected signs (fast)",{
  cond1 <- all(fast[prior1==0]==0)
  cond2 <- all(fast[prior1<0]<=0)
  cond3 <- all(fast[prior1>0]>=0)
  testthat::expect_true(cond1&cond2&cond3)
})

testthat::test_that("correlation (slow, fast)",{
  cond1 <- abs(mean(slow)-mean(fast))<0.01
  cond2 <- stats::cor(slow,fast)>0.99
  testthat::expect_true(cond1&cond2)
})

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#- - - equivalence predicted values- - - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for(scale in c("exp","iso")){
  
  family <- "gaussian"
  
  prior <- cbind(prior1,prior2)
  object <- transreg(y=y,X=X,prior=prior,family=family,scale=scale)
  y_hat1 <- predict(object,newx=X)
  
  beta <- coef(object$base,s=c(object$meta$lambda.min,object$meta$lambda.1se))
  omega <- as.numeric(coef(object$meta,s=object$meta$lambda.min))
  names <- paste0("source",seq_len(ncol(prior)))
  names(omega) <- c("(Intercept)",names,"lambda.min","lambda.1se")
  
  alpha_star <- omega["(Intercept)"] + omega["lambda.min"]*beta["(Intercept)","s1"] + omega["lambda.1se"]*beta["(Intercept)","s2"]
  beta_star <- rep(NA,times=p)
  
  if(scale=="exp"){
    for(j in seq_len(p)){
      beta_star[j] <- sum(omega[names]*object$base$prior$gamma*sign(object$base$z[j,])*abs(object$base$z[j,])^object$base$prior$tau) + omega["lambda.min"]*beta[1+j,"s1"] + omega["lambda.1se"]*beta[1+j,"s2"]
    }
  }
  if(scale=="iso"){
    for(j in seq_len(p)){
      beta_star[j] <- sum(omega[names]*object$base$prior$beta[j,]) + omega["lambda.min"]*beta[1+j,"s1"] + omega["lambda.1se"]*beta[1+j,"s2"]
    }
  }
  y_hat2 <- joinet:::.mean.function(alpha_star + X %*% beta_star,family=family)
  
  testthat::test_that("correlation (pred, coef)",{
    cond1 <- mean(y_hat1)-mean(y_hat2) < 0.01
    cond2 <- stats::cor(y_hat1,y_hat2) > 0.99
    testthat::expect_true(cond1&cond2)
  })
  
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#- - - prior re-scaling- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

family <- "gaussian"
prior <- object <- list()
prior$original <- prior1
prior$modified <- -stats::runif(1)*prior1

for(i in seq_along(prior)){
  set.seed(2)
  object[[i]] <- transreg(y=y,X=X,prior=prior[[i]],family=family,scale=scale)
}
pred <- list()
pred[[1]] <- predict(object[[1]],newx=X)
pred[[2]] <- predict(object[[2]],newx=X)

testthat::test_that("prior re-scaling",{
  cond1 <- mean(pred[[1]])-mean(pred[[2]])<0.01
  cond2 <- stats::cor(pred[[1]],pred[[2]])>0.99
  testthat::expect_true(cond1&cond2)
})

