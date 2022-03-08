
#--- slow versus fast isotonic scaling ---

n <- 100; p <- 1000
X <- matrix(stats::rnorm(n*p),nrow=n,ncol=p)
beta <- stats::rnorm(p)
prior <- beta + stats::rnorm(p)
y <- X %*% beta

slow <- iso.slow.single(y=y,X=X,prior=matrix(prior,ncol=1),family="gaussian")$beta
fast <- iso.fast.single(y=y,X=X,prior=matrix(prior,ncol=1),family="gaussian")$beta

table(prior=sign(prior),sign(slow)) # okay
table(prior=sign(prior),sign(fast)) # solve this problem
table(old=sign(slow),new=sign(fast))

testthat::test_that("expected signs (slow)",{
  cond1 <- all(slow[prior==0]==0)
  cond2 <- all(slow[prior<0]<=0)
  cond3 <- all(slow[prior>0]>=0)
  testthat::expect_true(cond1&cond2&cond3)
})

testthat::test_that("expected signs (fast)",{
  cond1 <- all(fast[prior==0]==0)
  cond2 <- all(fast[prior<0]<=0)
  cond3 <- all(fast[prior>0]>=0)
  testthat::expect_true(cond1&cond2&cond3)
})

testthat::test_that("correlation (slow, fast)",{
  cor <- stats::cor(slow,fast)
  testthat::expect_true(cor>0.99)
})


