library(svdcov)

source("data_helpers.R")

test_that("winsorize", {
  expect_equal(
    winsorize(-100:100, 0.9),
    c(rep(-90, 11), -89:89, rep(90, 11))
  )
})


test_that("winsorize matrix", {
  res <- winsorize(outer(-5:5, -5:5), 0.9)
  expected <- matrix(
    c(
      16, 16, 15, 10, 5, 0, -5, -10, -15, -16, -16,
      16, 16, 12, 8, 4, 0, -4, -8, -12, -16, -16,
      15, 12, 9, 6, 3, 0, -3, -6, -9, -12, -15,
      10, 8, 6, 4, 2, 0, -2, -4, -6, -8, -10,
      5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5,
      -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10,
      -15, -12, -9, -6, -3, 0, 3, 6, 9, 12, 15,
      -16, -16, -12, -8, -4, 0, 4, 8, 12, 16, 16,
      -16, -16, -15, -10, -5, 0, 5, 10, 15, 16, 16
    ),
    11, 11
  )
  expect_equal(res, expected)
})

test_that("winsorize error", {
  expect_error(
    winsorize(-100:100, 95)
  )
})

test_that("get_fitted", {
  svd_params <- SVDParams(
    row_cov=matrix(
      c(1, 0.5, 0.5,
        0.5, 1, 0.5,
        0.5, 0.5, 1),
      3, 3
    ),
    row_scores=data.frame(
      mean=c(1,2,3),
      score.1 = c(1, 1, -1),
      score.2 = c(-1, 1, 1)
    ),
    col_scores=data.frame(
      mean=c(1, 0, 1),
      score.1 = c(1, 1, -1),
      score.2 = c(-1, 1, 1)
    ),
    svd_d=c(2,1),
    log=FALSE
  )
  res <- get_fitted(svd_params)

  get_val <- function(row, col){
    svd_params@row_scores$mean[row] +
      svd_params@col_scores$mean[col] +
      sum(
        svd_params@svd_d *
          svd_params@row_scores[row, c("score.1", "score.2")] *
          svd_params@col_scores[col, c("score.1", "score.2")]
      )
  }

  for(i in 1:3) for(j in 1:3){
    expect_equal(res[i,j], get_val(i, j))
  }
})

zero_cov_data <- list(
  row_means=1:100,
  col_means=1:11,
  cov_mat=diag(100) * 0
)

zero_cov_data$x <- with(zero_cov_data, generate_cov_data(row_means, col_means, cov_mat))

test_that("Results on zero covariance, method=SVD", {
  results <- get_svd(zero_cov_data$x, method = "svd")

  real_means <- with(
    zero_cov_data,
    outer(row_means, col_means, FUN = `+`)
  )
  calculated_means <- get_fitted(results)
  expect_true(all(real_means == calculated_means))
})

test_that("Results on zero covariance, method=shrinkage", {
  results <- expect_warning(
    get_svd(zero_cov_data$x, method = "shrinkage")
  )
  real_means <- with(
    zero_cov_data,
    outer(row_means, col_means, FUN = `+`)
  )

  implied_means <- get_fitted(results)
  expect_true(all(real_means == implied_means))
})

test_that("block simulated_results, method=svd", {
  set.seed(215)
  data <- gen_block_diag_data(25, 0.3)
  results <- get_svd(data$x, method = "svd", winsorize=1.0, n_svd=4)

  real_means <- outer(data$row_means, data$col_means, FUN = `+`)
  implied_means <- get_fitted(results)
  if(FALSE){
    plot(real_means, implied_means)
    heatmap(data$cov_mat)
    heatmap(results@row_cov)
  }

  # This is actually small, since means are 1:100 :)
  expect_lt(mean(abs(real_means - implied_means)), 1)

  is_block <- get_mat_block(results@row_cov, 25)
  is_diag <- row(results@row_cov) == col(results@row_cov)

  diag_mean <- mean(diag(results@row_cov))
  block_mean <- mean(results@row_cov[is_block])
  offblock_mean <- mean(results@row_cov[!is_block & !is_diag])

  expect_equal(diag_mean - offblock_mean, 1.0, tolerance=0.05)
  expect_equal(block_mean - offblock_mean, 0.3, tolerance=0.05)
  expect_equal(diag_mean - block_mean, 0.7, tolerance=0.05)
})


test_that("block simulated_results, method=shrinkage", {
  set.seed(215)
  data <- gen_block_diag_data(25, 0.3)
  results <- get_svd(data$x, method = "shrinkage", winsorize=1.0, n_svd=4)

  real_means <- outer(data$row_means, data$col_means, FUN = `+`)
  implied_means <- get_fitted(results)
  if(FALSE){
    plot(real_means, implied_means)
    heatmap(data$cov_mat)
    heatmap(results@row_cov)
  }

  # 1 seems big, but is actually small, since means are 1:100 :)
  expect_lt(mean(abs(real_means - implied_means)), 1)

  is_block <- get_mat_block(results@row_cov, 25)
  is_diag <- row(results@row_cov) == col(results@row_cov)

  diag_mean <- mean(diag(results@row_cov))
  block_mean <- mean(results@row_cov[is_block])
  offblock_mean <- mean(results@row_cov[!is_block & !is_diag])

  expect_equal(diag_mean - offblock_mean, 1.0, tolerance=0.05)
  expect_equal(block_mean - offblock_mean, 0.3, tolerance=0.05)
  expect_equal(diag_mean - block_mean, 0.7, tolerance=0.05)
})


#
# test_that("sampled from posterior, block", {
#   set.seed(215)
#   method <- "shrinkage"
#   data <- block_diag_data$x[,1]
#   svd <- get_svd(block_diag_data$x[,-1], method = method, winsorize=FALSE, n_svd=4)
#   heatmap(svd@row_cov)
#
#   obs_id <- as.character(sample(seq_along(data), size = 10))
#   obs_data <- data[obs_id]
#
#   column_mean <- mean(obs_data - svd@row_scores$mean[as.numeric(obs_id)])
#
#   simulated_data <- sample_from_posterior(
#     svd,
#     obs_data,
#     obs_id,
#     column_mean=column_mean,
#     n_sim=400
#   )
#
#   sim_stats <- simulated_data %>%
#     filter(simulated) %>%
#     group_by(row) %>%
#     summarise(
#       mean = mean(value),
#       p025 = quantile(value, 0.025),
#       p975 = quantile(value, 0.975)
#     ) %>%
#     mutate(true_value = data[row])
#
#   plot(sim_stats$mean, sim_stats$true_value)
#   expect_close_to(
#     mean(sim_stats$mean - sim_stats$true_value),
#     0.0, 1.0
#   )
#
#   ci_covers <- (sim_stats$p025 < sim_stats$true_value) & (sim_stats$p975 > sim_stats$true_value)
#   expect_gt(mean(ci_covers), 0.8)
#
# })
#
# test_posterior_sample <- function(data_list, known_column_means, method){
#   data <- data_list$x[,1]
#   svd <- get_svd(
#     data_list$x[,-1],
#     method=method,
#     winsorize=FALSE,
#     known_column_means = known_column_means,
#     n_svd=4
#   )
#   heatmap(svd@row_cov)
#
#   obs_id <- as.character(sample(seq_along(data), size = 10))
#   obs_data <- data[obs_id]
#
#   if(is.null(known_column_means)){
#     column_mean <- mean(obs_data - svd@row_scores$mean[as.numeric(obs_id)])
#   } else{
#     column_mean <- known_column_means
#   }
#
#   simulated_data <- sample_from_posterior(
#     svd,
#     obs_data,
#     obs_id,
#     column_mean=column_mean,
#     n_sim=400
#   )
#
#   sim_stats <- simulated_data %>%
#     filter(simulated) %>%
#     group_by(row) %>%
#     summarise(
#       mean = mean(value),
#       p025 = quantile(value, 0.025),
#       p975 = quantile(value, 0.975)
#     ) %>%
#     mutate(true_value = data[row])
#
#   print(
#     ggplot(
#       sim_stats,
#       aes(x=mean, y=true_value)
#     ) + geom_point() +
#       geom_errorbarh(aes(xmin=p025, xmax=p975)) +
#       geom_abline(slope=1)
#   )
#
#   print(sprintf(
#     "Mean Deviation: %0.1f, Mean Abs: %0.1f, Mean True Value: %0.1f",
#     mean(sim_stats$mean - sim_stats$true_value),
#     mean(abs(sim_stats$mean - sim_stats$true_value)),
#     mean(sim_stats$true_value)
#   ))
#   expect_close_to(
#     mean(sim_stats$mean - sim_stats$true_value),
#     0.0, 1.0
#   )
#
#   print(sprintf("Cor(Pred, True): %0.2f", cor(sim_stats$mean, sim_stats$true_value)))
#   expect_gt(
#     cor(sim_stats$mean, sim_stats$true_value),
#     0.95
#   )
#
#   ci_covers <- (sim_stats$p025 < sim_stats$true_value) & (sim_stats$p975 > sim_stats$true_value)
#   print(sprintf("CI Covers: %s / %s = %0.2f", sum(ci_covers), length(ci_covers), mean(ci_covers)))
#   expect_gt(mean(ci_covers), 0.8)
# }
#
#
# test_that("sampled from posterior, block", {
#   set.seed(216)
#   test_posterior_sample(block_diag_data, known_column_means=NULL, method="svd")
# })
#
#
# test_that("sampled from posterior, block", {
#   set.seed(217)
#   test_posterior_sample(block_diag_data, known_column_means=0, method="svd")
# })
#
# test_that("sampled from posterior, block",{
#   set.seed(218)
#   test_posterior_sample(block_diag_data, known_column_means=NULL, method="shrinkage")
# })
#
# test_that("sampled from posterior, block",{
#   set.seed(219)
#   test_posterior_sample(block_diag_data, known_column_means=0, method="shrinkage")
# })
