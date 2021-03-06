
# svdcov

<!-- badges: start -->
<!-- badges: end -->

The goal of `svdcov` is to calculate the covariance matrix of rectangular data 
using SVD or shrinkage methods.

## Installation

You can install the released version of svdcov from github with:

``` r
devtools::intall_github("github.io/jtannen/svdcov")
```

## Example


``` r
library(svdcov)

# Consider a block-diagonal example.
cov <- matrix(0, 10, 10)
cov[row(cov) <=5 & col(cov) <= 5] <- 0.5
cov[row(cov) > 5 & col(cov) > 5] <- 0.5
diag(cov) <- 1

data <- MASS::mvrnorm(
  n=100, 
  mu=1:10, 
  Sigma=cov
)

data <- t(data)

res <- get_svd(data, n_svd=2, col_means=rep(0, 100))

res@row_cov
get_fitted(res)

sample_from_posterior(
  res,
  obs=c(1.42, 7.34),
  obs_id=c(1, 6),
  col_mean=0,
  col_mean_sd=1
)
```

