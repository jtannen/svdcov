generate_cov_data <- function(
  row_means,
  col_means,
  cov_mat
){
  x <- t(
    MASS::mvrnorm(
      n = length(col_means),
      mu = row_means,
      Sigma = cov_mat
    ) + col_means
  )

  row.names(x) <- seq_along(row_means)
  colnames(x) <- seq_along(col_means)
  x
}


get_mat_diag <- function(mat){
  return(row(mat) == col(mat))
}

get_mat_block <- function(mat, blocksize){
  row_mod <- (row(mat)-1) %/% blocksize
  col_mod <- (col(mat)-1) %/% blocksize
  is_diag <- get_mat_diag(mat)

  is_block <- (row_mod == col_mod) & !is_diag
  return(is_block)
}

generate_block_cov_mat <- function(nrow = 100, blocksize=25, offdiag=0.9){
  block_cov <- diag(nrow)
  is_block <- get_mat_block(block_cov, blocksize)
  block_cov[is_block] <- offdiag
  block_cov
}

gen_block_diag_data <- function(blocksize, offdiag){
  block_diag_data <- list(
    row_means = 1:100,
    col_means = 1:200
  )
  block_diag_data$cov_mat <-
    generate_block_cov_mat(length(block_diag_data$row_means), blocksize=blocksize, offdiag=offdiag)

  block_diag_data$x <- with(
    block_diag_data,
    generate_cov_data(row_means, col_means, cov_mat)
  )
  block_diag_data
}

