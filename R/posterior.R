
mvrnorm_chol <- function(n, mu, sigma_chol){
  sample_raw <- matrix(rnorm(n * length(mu)), length(mu), n)
  sample <- t(sigma_chol) %*% sample_raw + matrix(mu, nrow=length(mu), ncol=n)
  rownames(sample) <- rownames(sigma_chol)
  return(t(sample))
}

mvrnorm_chol_mat <- function(mu_mat, sigma_chol){
  sample_raw <- matrix(
    rnorm(prod(dim(mu_mat))),
    nrow(mu_mat),
    ncol(mu_mat)
  )
  sample <- t(sigma_chol) %*% sample_raw + mu_mat
  rownames(sample) <- rownames(sigma_chol)
  return(t(sample))
}

memoised_chol <- memoise::memoise(chol)


#' Sample a new column from an SVDParams, using partially observed row values.
#'
#' @description
#' Suppose D is a NxM matrix of historical observations.
#' `svd` is the SVDParams of the historical data.
#' This samples a new column for the data with covariance svd@row_cov.
#' You can also specify (a) a mean and sd for the col_mean,
#' and (b) a subset of the values that have been observed (obs and obs_id).
#'

sample_from_posterior <- function(
  svd,
  obs,
  obs_id,
  col_mean_prior=NULL,
  col_mean_prior_var=NULL,
  filter_to_ids=NULL,
  n_sim=400,
  verbose=FALSE
){
  vprint <- function(x) if(verbose) print(x) else NULL

  row_cov <- svd@row_cov
  row_scores <- svd@row_scores

  if(!is.null(filter_to_ids)){
    use_rows <- row.names(row_cov) %in% filter_to_ids
    row_cov <- row_cov[use_rows, use_rows]
    row_scores <- row_scores %>% filter(row %in% filter_to_ids)
    use_obs <- obs_id %in% filter_to_ids
    obs <- obs[use_obs]
    obs_id <- obs_id[use_obs]
  }

  row_means <- row_scores$mean
  row_ids <- row_scores$row

  vprint("Calculating Matrices")

  a_rows <- match(obs_id, row_ids)
  c_rows <- seq_along(row_means)[-a_rows]

  ## Assume distribution is normal with covariance [A, B], [B', C]
  A <- row_cov[a_rows, a_rows] # N_A x N_A
  B <- row_cov[c_rows, a_rows] # N_C x N_A
  C <- row_cov[c_rows, c_rows] # N_C x N_C

  a <- obs # length N_A

  vprint("Calculating A inverse")
  vprint(sprintf("nrow(A) = %s", nrow(A)))

  a_dm <- a - row_means[a_rows]
  if(is.null(col_mean_prior)){
    col_mean_prior <- mean(svd@col_scores$mean)
    col_mean_prior_var <- var(svd@col_scores$mean)
  }

  col_pars <- posterior_mean(a_dm, A, col_mean_prior, col_mean_prior_var)

  sim_means <- rnorm(n_sim, mean=col_pars$mean, sd=sqrt(col_pars$var)) # n_sim
  prior_means <- outer(row_means, sim_means, FUN="+") # (N_A + N_C) x n_sim

  if(length(a_rows) > 0){
    A_inv <- solve(A)
  } else {
    A_inv <- A
  }

  c_means <- (
    prior_means[c_rows,] # N_C x n_sim
    + B %*% A_inv %*% (a - prior_means[a_rows,]) # N_C x n_sim
  )

  c_sigma <- C - B %*% A_inv %*% t(B) # N_C x N_C
  if(length(a_rows) <= 1) rownames(c_sigma) <- names(B) ## otherwise it gets inherited

  vprint("Calculating Cholesky")

  vprint("MVRNorm")
  if(nrow(c_sigma) > 0){
    chol_c_sigma <- memoised_chol(c_sigma) # N_C x N_C
    sim <- mvrnorm_chol_mat(c_means, sigma_chol=chol_c_sigma)
    sim <- as.data.frame(sim) %>%
      mutate(sim = 1:n_sim) %>%
      gather("row", "value", -sim)
  } else {
    sim <- data.frame(sim = integer(0), row = character(0), value = numeric(0))
  }

  vprint("Processing Simulations")

  obs_rep <- bind_rows(
    rep(
      list(data.frame(row=as.character(obs_id), value=obs)),
      n_sim
    ),
    .id="sim"
  ) %>% mutate(simulated=FALSE)

  sim <- sim %>%
    mutate(simulated=TRUE, sim=as.character(sim)) %>%
    bind_rows(obs_rep)

  return(sim)
}

posterior_mean <- function(x, Sigma, mu_0, var_0){
  # suppose
  # x ~ norm(mu, Sigma)
  # mu ~ norm(mu_0, var_0)
  # Then LL(x | mu) = - 0.5 (x - 1mu)' Sigma^-1 (x - 1mu) - 0.5 (mu - mu_0)^2/var_0
  # So
  # var_post = (Sigma^-1 + 1/var_0)^-1
  # mu | x ~ norm(var_post * (Sigma^-1 x + mu_0 / var_0), var_post)

  var_post <- (sum(solve(Sigma)) + 1/var_0)^-1
  mean_post <- var_post * (sum(solve(Sigma, x)) + mu_0 / var_0)

  list(
    mean=mean_post,
    var=var_post
  )
}


posterior_sample <- function(...) stop("Deprecated. Change to sample_from_posterior().")


sample_from_posterior_backup <- function(
  svd,
  obs,
  obs_id,
  column_mean,  ## prior for column_mean
  filter_to_ids=NULL,
  n_sim=400,
  verbose=FALSE
){

  vprint <- function(x) if(verbose) print(x) else NULL

  row_cov <- svd@row_cov
  row_scores <- svd@row_scores

  if(!is.null(filter_to_ids)){
    use_rows <- row.names(row_cov) %in% filter_to_ids
    row_cov <- row_cov[use_rows, use_rows]
    row_scores <- row_scores %>% filter(row %in% filter_to_ids)
    use_obs <- obs_id %in% filter_to_ids
    obs <- obs[use_obs]
    obs_id <- obs_id[use_obs]
  }

  row_means <- row_scores$mean
  row_ids <- row_scores$row

  vprint("Calculating Matrices")

  a_rows <- match(obs_id, row_ids)
  c_rows <- seq_along(row_means)[-a_rows]

  ## Assume distribution is normal with covariance [A, B], [B', C]
  A <- row_cov[a_rows, a_rows] # N_A x N_A
  B <- row_cov[c_rows, a_rows] # N_C x N_A
  C <- row_cov[c_rows, c_rows] # N_C x N_C

  a <- obs # N_A

  vprint("Calculating A inverse")
  vprint(sprintf("nrow(A) = %s", nrow(A)))

  if(length(a_rows) > 0){
    A_inv <- solve(A)
  } else {
    A_inv <- A
  }

  prior_mean <- row_means + column_mean

  c_mean <- prior_mean[c_rows] + B %*% A_inv %*% (a - prior_mean[a_rows])

  c_sigma <- C - B %*% A_inv %*% t(B)
  if(length(a_rows) <= 1) rownames(c_sigma) <- names(B) ## otherwise it gets inherited

  vprint("Calculating Cholesky")
  chol_c_sigma <- memoised_chol(c_sigma)

  vprint("MVRNorm")
  if(nrow(c_sigma) > 0){
    sim <- mvrnorm_chol(n=n_sim, mu=c_mean, sigma_chol=chol_c_sigma)
    sim <- as.data.frame(sim) %>%
      mutate(sim = 1:n_sim) %>%
      gather("row", "value", -sim)
  } else {
    sim <- data.frame(sim = integer(0), row = character(0), value = numeric(0))
  }

  vprint("Processing Simulations")

  obs_rep <- bind_rows(
    rep(
      list(data.frame(row=obs_id, value=obs)),
      n_sim
    ),
    .id="sim"
  ) %>% mutate(simulated=FALSE)

  sim <- sim %>%
    mutate(simulated=TRUE, sim=as.character(sim)) %>%
    bind_rows(obs_rep)

  return(sim)
}

