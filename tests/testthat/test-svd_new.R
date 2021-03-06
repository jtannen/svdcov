# library(svdcov)
#
# test_that("winsorize", {
#   expect_equal(
#     winsorize(-100:100, 0.9),
#     c(rep(-90, 11), -89:89, rep(90, 11))
#   )
# })
#
#
# test_that("winsorize matrix", {
#   res <- winsorize(outer(-5:5, -5:5), 0.9)
#   expected <- matrix(
#     c(
#       16, 16, 15, 10, 5, 0, -5, -10, -15, -16, -16,
#       16, 16, 12, 8, 4, 0, -4, -8, -12, -16, -16,
#       15, 12, 9, 6, 3, 0, -3, -6, -9, -12, -15,
#       10, 8, 6, 4, 2, 0, -2, -4, -6, -8, -10,
#       5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5,
#       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#       -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5,
#       -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10,
#       -15, -12, -9, -6, -3, 0, 3, 6, 9, 12, 15,
#       -16, -16, -12, -8, -4, 0, 4, 8, 12, 16, 16,
#       -16, -16, -15, -10, -5, 0, 5, 10, 15, 16, 16
#     ),
#     11, 11
#   )
#   expect_equal(res, expected)
# })
#
# test_that("winsorize error", {
#   expect_error(
#     winsorize(-100:100, 95)
#   )
# })
#
# test_that("get_fitted", {
#   svd_params <- SVDParams(
#     row_cov <- matrix(
#       c(1, 0.5, 0.5,
#         0.5, 1, 0.5,
#         0.5, 0.5, 1),
#       3, 3
#     ),
#     row_scores=data.frame(
#       mean=c(1,2,3),
#       score.1 = c(1, 1, -1),
#       score.2 = c(-1, 1, 1)
#     ),
#     col_scores=data.frame(
#       mean=c(1, 0, 1),
#       score.1 = c(1, 1, -1),
#       score.2 = c(-1, 1, 1)
#     ),
#     svd_d=c(2,1),
#     log=FALSE
#   )
#   res <- get_fitted(svd_params)
#
#   get_val <- function(row, col){
#     svd_params@row_scores$mean[row] +
#       svd_params@col_scores$mean[col] +
#       sum(
#         svd_params@svd_d *
#           svd_params@row_scores[row, c("score.1", "score.2")] *
#           svd_params@col_scores[col, c("score.1", "score.2")]
#       )
#   }
#
#   for(i in 1:3) for(j in 1:3){
#     expect_equal(res[i,j], get_val(i, j))
#   }
# })
#
#
# generate_cov_data <- function(
#   row_means,
#   col_means,
#   cov_mat
# ){
#   x <- t(
#     MASS::mvrnorm(
#       n = length(col_means),
#       mu = row_means,
#       Sigma = cov_mat
#     ) + col_means
#   )
#
#   row.names(x) <- seq_along(row_means)
#   colnames(x) <- seq_along(col_means)
#   x
# }
#
# zero_cov_data <- list(
#   row_means=1:100,
#   col_means=1:11,
#   cov_mat=diag(100) * 0
# )
#
# zero_cov_data$x <- with(zero_cov_data, generate_cov_data(row_means, col_means, cov_mat))
#
# get_mat_diag <- function(mat){
#   return(row(mat) == col(mat))
# }
#
# get_mat_block <- function(mat, blocksize){
#   row_mod <- (row(mat)-1) %/% blocksize
#   col_mod <- (col(mat)-1) %/% blocksize
#   is_diag <- get_mat_diag(mat)
#
#   is_block <- (row_mod == col_mod) & !is_diag
#   return(is_block)
# }
#
# generate_block_cov_mat <- function(nrow = 100, blocksize=25, offdiag=0.9){
#   block_cov <- diag(nrow)
#   is_block <- get_mat_block(block_cov, blocksize)
#   block_cov[is_block] <- offdiag
#   block_cov
# }
#
# extract_mean_mat <- function(svd_results){
#   outer(
#     svd_results@row_scores$mean,
#     svd_results@col_scores$mean,
#     FUN = `+`
#   )
# }
#
# test_that("Results on zero covariance, method=SVD", {
#   results <- get_svd(zero_cov_data$x, method = "svd")
#
#   real_means <- with(
#     zero_cov_data,
#     outer(row_means, col_means, FUN = `+`)
#   )
#   calculated_means <- get_fitted(results)
#   expect_true(all(real_means == calculated_means))
# })
#
# test_that("Results on zero covariance, method=shrinkage", {
#   results <- expect_warning(
#     get_svd(zero_cov_data$x, method = "shrinkage")
#   )
#   real_means <- with(
#     zero_cov_data,
#     outer(row_means, col_means, FUN = `+`)
#   )
#
#   implied_means <- get_fitted(results)
#   expect_true(all(real_means == implied_means))
# })
#
# gen_block_diag_data <- function(blocksize, offdiag){
#   block_diag_data <- list(
#     row_means = 1:100,
#     col_means = 1:200
#   )
#   block_diag_data$cov_mat <-
#     generate_block_cov_mat(length(block_diag_data$row_means), blocksize=blocksize, offdiag=offdiag)
#
#   block_diag_data$x <- with(
#     block_diag_data,
#     generate_cov_data(row_means, col_means, cov_mat)
#   )
#   block_diag_data
# }
#
#
# test_that("block simulated_results, method=svd", {
#   set.seed(215)
#   data <- gen_block_diag_data(25, 0.3)
#   results <- get_svd(data$x, method = "svd", winsorize=1.0, n_svd=4)
#
#   real_means <- outer(data$row_means, data$col_means, FUN = `+`)
#   implied_means <- get_fitted(results)
#   if(FALSE){
#     plot(real_means, implied_means)
#     heatmap(data$cov_mat)
#     heatmap(results@row_cov)
#   }
#
#   # This is really small, since means are 1:100 :)
#   expect_lt(mean(abs(real_means - implied_means)), 1)
#
#   is_block <- get_mat_block(results@row_cov, 25)
#   is_diag <- row(results@row_cov) == col(results@row_cov)
#
#   diag_mean <- mean(diag(results@row_cov))
#   block_mean <- mean(results@row_cov[is_block])
#   offblock_mean <- mean(results@row_cov[!is_block & !is_diag])
#
#   expect_equal(diag_mean - offblock_mean, 1.0, tolerance=0.05)
#   expect_equal(block_mean - offblock_mean, 0.3, tolerance=0.05)
#   expect_equal(diag_mean - block_mean, 0.7, tolerance=0.05)
# })
#
#
# test_that("block simulated_results, method=shrinkage", {
#   set.seed(215)
#   data <- gen_block_diag_data(25, 0.3)
#   results <- get_svd(data$x, method = "shrinkage", winsorize=1.0, n_svd=4)
#
#   real_means <- outer(data$row_means, data$col_means, FUN = `+`)
#   implied_means <- get_fitted(results)
#   if(FALSE){
#     plot(real_means, implied_means)
#     heatmap(data$cov_mat)
#     heatmap(results@row_cov)
#   }
#
#   # 1 seems big, but is actually small, since means are 1:100 :)
#   expect_lt(mean(abs(real_means - implied_means)), 1)
#
#   is_block <- get_mat_block(results@row_cov, 25)
#   is_diag <- row(results@row_cov) == col(results@row_cov)
#
#   diag_mean <- mean(diag(results@row_cov))
#   block_mean <- mean(results@row_cov[is_block])
#   offblock_mean <- mean(results@row_cov[!is_block & !is_diag])
#
#   expect_equal(diag_mean - offblock_mean, 1.0, tolerance=0.05)
#   expect_equal(block_mean - offblock_mean, 0.3, tolerance=0.05)
#   expect_equal(diag_mean - block_mean, 0.7, tolerance=0.05)
# })



