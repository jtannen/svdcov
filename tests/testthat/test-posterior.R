library(svdcov)

source("data_helpers.R")

test_that("mvrnorm_chol", {
  set.seed(215)
  cov <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  chol_cov <- chol(cov)

  res <- mvrnorm_chol(1e6, 1:2, chol_cov)
  expect_equal(dim(res), c(1e6, 2))
  expect_equal(colMeans(res), 1:2, tolerance = 0.005)
  expect_equal(cov(res), cov, tolerance = 0.005)
})

test_that("mvrnorm_chol_mat", {
  set.seed(215)
  cov <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  chol_cov <- chol(cov)

  means <- matrix(1:2e6, nrow=2, byrow = T)

  res <- mvrnorm_chol_mat(means, chol_cov)
  expect_equal(dim(res), c(1e6, 2))
  expect_equal(cor(means[1,], res[,1]), 1, tolerance=1e-6)
  expect_equal(cor(means[2,], res[,2]), 1, tolerance=1e-6)
  expect_equal(colMeans(res), rowMeans(means), tolerance = 1e-6)
  expect_equal(cov(res - t(means)), cov, tolerance = 0.005)
})


test_that("posterior_mean", {
  Sigma <- matrix(0, 1, 1)
  diag(Sigma) <- 1
  mu_0 <- 0
  var_0 <- 1
  x <- c(1)
  res <- posterior_mean(x, Sigma, mu_0, var_0)
  expect_equal(res$mean, 0.5)
  expect_equal(res$var, 0.5)

  Sigma <- matrix(0, 2, 2)
  diag(Sigma) <- 1
  mu_0 <- 0
  var_0 <- 1
  x <- c(1, 1)
  res <- posterior_mean(x, Sigma, mu_0, var_0)
  expect_equal(res$mean, 2/3)
  expect_equal(res$var, 1/3)

  Sigma <- matrix(0, 2, 2)
  diag(Sigma) <- 1
  mu_0 <- 1
  var_0 <- 1
  x <- c(1, 1)
  res <- posterior_mean(x, Sigma, mu_0, var_0)
  expect_equal(res$mean, 1)
  expect_equal(res$var, 1/3)

  Sigma <- matrix(0.5, 2, 2)
  diag(Sigma) <- 1
  mu_0 <- 0
  var_0 <- 1
  x <- c(1, 1)
  res <- posterior_mean(x, Sigma, mu_0, var_0)
  expect_equal(res$mean, 4/7)
  expect_equal(res$var, 3/7)

})

test_that("block sample_from_posterior, null column prior", {
  set.seed(215)

  # No scores since those are ignored
  svd_params <- SVDParams(
    row_cov = matrix(
      c(1, 0.3, 0, 0, 0,
        0.3, 1, 0, 0, 0,
        0, 0, 1, 0.3, 0.3,
        0, 0, 0.3, 1, 0.3,
        0, 0, 0.3, 0.3, 1),
      5, 5
    ),
    row_scores=data.frame(
      row=1:5,
      mean=1:5
    ),
    col_scores=data.frame(
      column=1:3,
      mean=c(10, 11, 9) #mean=10, var=1
    ),
    svd_d=numeric(),
    log=FALSE
  )

  res <- sample_from_posterior(
    svd=svd_params,
    obs=c(10, 14), # -1 for row 1, +1 for row 3
    obs_id=as.character(c(1, 3)),
    n_sim = 10e3
  )

  res_summ <- res %>% group_by(row) %>%
    summarise(
      is_sim = mean(simulated),
      value_mean = mean(value),
      value_var = var(value),
      .groups="drop"
    )
  expect_equal(
    res_summ$is_sim,
    c(0, 1, 0, 1, 1)
  )

  expect_equal(
    res_summ$value_mean,
    10 + 1:5 + c(-1, -0.3, 1, 0.3, 0.3),
    tolerance = 0.05
  )

  res_mat <- matrix(NA, max(as.numeric(res$sim)), max(as.numeric(res$row)))
  res_mat[cbind(as.numeric(res$sim), as.numeric(res$row))] <- res$value
  cov_res <- cov(res_mat)

  var_from_col_means <- 0.333 * (1-0.3)^2
  var_from_chol <- (1 - 0.3^2)
  cov_from_chol <- (0.3 - 0.3^2)
  total_var <- var_from_col_means + var_from_chol

  expect_equal(
    diag(cov_res),
    c(0, total_var, 0, total_var, total_var),
    tolerance = 0.02
  )

  expect_equal(
    cov_res[2, 4:5],
    rep(var_from_col_means, 2),
    tolerance = 0.02
  )

  expect_equal(
    cov_res[4, 5],
    var_from_col_means + cov_from_chol,
    tolerance = 0.05
  )

})


test_that("block sample_from_posterior, known column prior", {
  set.seed(215)

  # No scores since those are ignored
  svd_params <- SVDParams(
    row_cov = matrix(
      c(1, 0.3, 0, 0, 0,
        0.3, 1, 0, 0, 0,
        0, 0, 1, 0.3, 0.3,
        0, 0, 0.3, 1, 0.3,
        0, 0, 0.3, 0.3, 1),
      5, 5
    ),
    row_scores=data.frame(
      row=1:5,
      mean=1:5
    ),
    col_scores=data.frame(
      column=1:3,
      mean=c(10, 11, 9) #mean=10, var=1
    ),
    svd_d=numeric(),
    log=FALSE
  )

  res <- sample_from_posterior(
    svd=svd_params,
    obs=c(10, 14), # -1 for row 1, +1 for row 3
    obs_id=as.character(c(1, 3)),
    col_mean_prior=11,
    col_mean_prior_var=2,
    n_sim = 10e3
  )

  res_summ <- res %>% group_by(row) %>%
    summarise(
      is_sim = mean(simulated),
      value_mean = mean(value),
      value_var = var(value),
      .groups="drop"
    )
  expect_equal(
    res_summ$is_sim,
    c(0, 1, 0, 1, 1)
  )

  expect_equal(
    res_summ$value_mean,
    10 + 1:5 + c(-1, -0.3, 1, 0.3, 0.3),
    tolerance = 0.05
  )

  res_mat <- matrix(NA, max(as.numeric(res$sim)), max(as.numeric(res$row)))
  res_mat[cbind(as.numeric(res$sim), as.numeric(res$row))] <- res$value
  cov_res <- cov(res_mat)

  var_from_col_means <- 0.4 * (1-0.3)^2
  var_from_chol <- (1 - 0.3^2)
  cov_from_chol <- (0.3 - 0.3^2)
  total_var <- var_from_col_means + var_from_chol

  expect_equal(
    diag(cov_res),
    c(0, total_var, 0, total_var, total_var),
    tolerance = 0.02
  )

  expect_equal(
    cov_res[2, 4:5],
    rep(var_from_col_means, 2),
    tolerance = 0.02
  )

  expect_equal(
    cov_res[4, 5],
    var_from_col_means + cov_from_chol,
    tolerance = 0.05
  )

})
