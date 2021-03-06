#' @import dplyr

get_svd <- function(
  mat,
  n_svd=3,
  winsorize=1.0, # set to lower to winsorize
  row_groups=NULL,
  col_means=NULL,  ## Priors for mean. If not null, use fixed effects.
  verbose=FALSE,
  method=c("svd", "shrinkage"),
  is_logged=FALSE
){
  method=match.arg(method)

  if(is.null(row.names(mat))) row.names(mat) <- 1:nrow(mat)
  if(is.null(colnames(mat))) colnames(mat) <- 1:ncol(mat)

  if(is.null(row_groups)){
    row_groups <- 1:nrow(mat)
  } else {
    assertthat::assert_that(length(row_groups) == nrow(mat))
  }

  row_means <- rowMeans(mat)
  mat_demeaned <- sweep(mat, row_means, MARGIN=1)

  if(is.null(col_means)){
    col_means <- colMeans(mat_demeaned)
  } else {
    assertthat::assert_that(length(col_means) == ncol(mat))
  }

  mat_demeaned <- sweep(mat_demeaned, col_means, MARGIN=2)

  mat_demeaned <- winsorize(mat_demeaned, prob=winsorize)

  if(method == "shrinkage"){
    svd_params <- get_shrinkage_cov(mat_demeaned, row_means, col_means, n_svd, verbose=verbose)
  } else if (method == "svd"){
    svd_params <- get_svd_cov(mat_demeaned, row_means, col_means, n_svd, row_groups)
  } else {
    stop("method must be ('shrinkage', 'svd')")
  }

  svd_params@log <- is_logged

  if(verbose){
    stop("Not implemented.")
    # print_fitted_vs_obs(svd_params, mat_demeaned)
    # print_ggplot_scatter(svd_params@row_cov, mat_demeaned)
  }

  return(svd_params)
}

winsorize <- function(x, prob){
  if(prob <= 0) stop("Bad winsorization prob. Did you set as FALSE? Use 1.0 instead.")
  if(prob > 1) stop("winsorization prob needs to be in (0, 1).")

  threshold <- quantile(abs(x), prob)
  replace <- abs(x) > threshold
  x[replace] <- sign(x[replace]) * threshold
  x
}

get_fitted <- function(svd){
  row_scores <- svd@row_scores
  col_scores <- svd@col_scores
  svd_d <- svd@svd_d
  row_score_mat <- as.matrix(row_scores %>% select(starts_with("score.")))
  col_score_mat <- as.matrix(col_scores %>% select(starts_with("score.")))

  means <- outer(row_scores$mean, col_scores$mean, FUN = `+`)
  fitted <- means + row_score_mat %*% diag(svd_d) %*% t(col_score_mat)

  if(svd@log) fitted <- exp(fitted)

  return(fitted)
}


## hacky way to convert cov.shrink type to matrix
shrinkage_cov_as_matrix <- function(shrinkage){
  k <- nrow(shrinkage)
  if(k != ncol(shrinkage)) stop("shrinkage isn't square. It should be.")
  return(shrinkage[1:k, 1:k])
}

get_shrinkage_cov <- function(
  mat_demeaned,
  row_means,
  col_means,
  n_svd,
  verbose,
  ...
){
  ## Uses corpcor::cov.shrink to calculate historic covariance

  row_cov <- corpcor::cov.shrink(t(mat_demeaned), verbose=verbose)
  row_cov <- shrinkage_cov_as_matrix(row_cov)
  e <- eigen(row_cov, symmetric=TRUE)

  svd_d <- e$values[1:n_svd]

  row_scores <- data.frame(
    row = row.names(mat_demeaned),
    score = e$vectors[,1:n_svd],
    mean = row_means
  )

  ## col scores are dot product between row_scores and obs
  col_score_mat <- t(mat_demeaned) %*% e$vectors[, 1:n_svd] %*% diag(1/svd_d)
  col_score_mat[ ,svd_d == 0] <- 0

  col_scores <- data.frame(
    col = colnames(mat_demeaned),
    score = col_score_mat,
    mean = col_means
  )

  return(
    SVDParams(
      row_cov=row_cov,
      row_scores=row_scores,
      col_scores=col_scores,
      svd_d=svd_d
    )
  )
}


get_svd_cov <- function(
  mat_demeaned,
  row_means,
  col_means,
  n_svd,
  row_groups,
  verbose=TRUE,
  ...
){
  ## Uses the leading svd dimensions to calculate historic covariance
  ## This shrinks much less than cov.shrink

  svd <- svd(mat_demeaned, n_svd, n_svd)
  svd_d <- svd$d[1:n_svd]

  fitted <- svd$u %*% diag(svd_d) %*% t(svd$v)
  resid <- (mat_demeaned - fitted) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("row") %>%
    gather("column", "residual", -row)

  resid_var <- var(as.vector(fitted - mat_demeaned))

  vprint <- function(...) if(verbose) print(...)
  vprint("Total Resid Var:"); vprint(resid_var)

  row_cov <- svd$u %*% diag(svd_d) %*% cov(svd$v) %*% diag(svd_d) %*% t(svd$u)
  diag(row_cov) <- diag(row_cov) + resid_var

  if(!all(diff(row_groups)==1)){
    vprint("Calculating group variance")
    group_resid <- resid %>%
      mutate(group=row_groups) %>%
      group_by(group, column) %>%
      summarise(mean_resid = mean(residual))

    group_var <- var(group_resid$mean_resid)
    vprint("Group Variance:"); vprint(group_var)

    for(group in unique(row_groups)){
      rows <- which(row_groups == group)
      row_cov[rows, rows] <- row_cov[rows, rows] + group_var
    }
    diag(row_cov) <- diag(row_cov) - group_var
  }

  row_scores <- data.frame(
    row = row.names(mat_demeaned),
    score = svd$u,
    mean = row_means
  )

  col_scores <- data.frame(
    col = colnames(mat_demeaned),
    score = svd$v,
    mean = col_means
  )

  return(
    SVDParams(
      row_cov=row_cov,
      row_scores=row_scores,
      col_scores=col_scores,
      svd_d=svd_d
    )
  )
}
