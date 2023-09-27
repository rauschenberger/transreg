
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

fast <- .iso.fast.single(y=y,X=X,prior=matrix(prior1,ncol=1),family="gaussian")$beta
  
testthat::test_that("expected signs (fast)",{
  cond1 <- all(fast[prior1==0]==0)
  cond2 <- all(fast[prior1<0]<=0)
  cond3 <- all(fast[prior1>0]>=0)
  testthat::expect_true(cond1&cond2&cond3)
})

if(require("CVXR")){
  
  slow <- .iso.slow.single(y=y,X=X,prior=matrix(prior1,ncol=1),family="gaussian")$beta

  testthat::test_that("expected signs (slow)",{
    cond1 <- all(slow[prior1==0]==0)
    cond2 <- all(slow[prior1<0]<=0)
    cond3 <- all(slow[prior1>0]>=0)
    testthat::expect_true(cond1&cond2&cond3)
  })
  
  testthat::test_that("correlation (slow, fast)",{
    cond1 <- abs(mean(slow)-mean(fast))<0.01
    cond2 <- stats::cor(slow,fast)>0.99
    testthat::expect_true(cond1&cond2)
  })
  
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#- - - equivalence predicted values- - - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for(scale in c("exp","iso")){
  for(stack in c("sta","sim")){
    
    family <- "gaussian"
    
    prior <- cbind(prior1,prior2)
    object <- transreg(y=y,X=X,prior=prior,family=family,scale=scale,stack=stack)
    
    y_hat1 <- predict(object,newx=X)
    coef <- coef(object=object)
    y_hat2 <- joinet:::.mean.function(coef$alpha + X %*% coef$beta,family=family)
    
    testthat::test_that("correlation (pred, coef)",{
      cond1 <- mean(y_hat1)-mean(y_hat2) < 0.01
      cond2 <- stats::cor(y_hat1,y_hat2) > 0.99
      testthat::expect_true(cond1&cond2)
    })
  
  }
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#- - - prior re-scaling- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

family <- "gaussian"
prior <- object <- pred <- list()
prior$original <- prior1
prior$modified <- -stats::runif(1)*prior1

for(scale in c("exp","iso")){
  for(stack in c("sta","sim")){
    
    for(i in seq_along(prior)){
      set.seed(2)
      object[[i]] <- transreg(y=y,X=X,prior=prior[[i]],family=family,scale=scale,stack=stack,switch=TRUE)
      pred[[i]] <- predict(object[[i]],newx=X)
    }
    
    testthat::test_that("prior re-scaling",{
      cond1 <- mean(pred[[1]])-mean(pred[[2]])<0.01
      cond2 <- stats::cor(pred[[1]],pred[[2]])>0.99
      testthat::expect_true(cond1&cond2)
    })
    
  }
}

