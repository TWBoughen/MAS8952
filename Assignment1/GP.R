# Samples from MVN
rmvn_svd = function(N, mean, covariance) {
  # samples from the MVN using singular value decomposition
  
  n <- length(mean) # dimension of random vector
  
  standard_normal_sample <- rnorm(n * N)
  standard_mvn_sample <- matrix(standard_normal_sample, nrow = N)
  
  svd_out <- svd(covariance)
  square_root =  svd_out$v %*% sqrt(diag(svd_out$d))
  output <- square_root %*% t(standard_mvn_sample) + mean
  
  return(t(output))
}

# Computes the log density of MVN
dmvn = function(x, mean, covariance) {
  # computes log density of MVN
  if (is.vector(x)) {
    x = matrix(x, ncol = length(x))
  }
  n = ncol(x)
  
  cholesky = chol(covariance)
  solved_cholesky = backsolve(cholesky, t(x) - mean, transpose = TRUE)
  quad_form = colSums(solved_cholesky**2)
  output = -sum(log(diag(cholesky))) - 0.5 * n * log(2 * pi) - 0.5 * quad_form
  names(output) = rownames(x)
  return(output)
}

# Computes the conditional MVN
condmvn = function(a, mu1, mu2, K11, K12, K22) {
  # computes conditional of MVN
  
  inv_K22 = solve(K22)
  
  tilde_mu = mu1 + K12 %*% inv_K22 %*% (a - mu2)
  tilde_K = K11 - K12 %*% inv_K22 %*% t(K12)
  
  output = list("mean_vec" = as.vector(tilde_mu), "cov_mat" = tilde_K)
  
  return(output)
}

# Computes posterior Gaussian Process
gp_posterior = function(x_pred, x_cond, y_cond, mean_func, kernel_func, error_var = 0, ...) {
  # Computes the posterior mean vector and posterior covariance matrix
  # at locations x_pred from
  # conditioning GP f = y_cond at locations x_cond
  
  cov_mat_cond = kernel_func(x_cond, x_cond, ...)
  cov_mat_pred = kernel_func(x_pred, x_pred, ...)
  cross_covariance = kernel_func(x_pred, x_cond, ...)
  
  n = nrow(cov_mat_cond)
  error_cov_mat = diag(n) * error_var
  cov_mat_cond = cov_mat_cond + error_cov_mat
  
  mean_cond = mean_func(x_cond, ...)
  mean_pred = mean_func(x_pred, ...)
  
  condmvn_output = condmvn(y_cond, mean_pred, mean_cond, cov_mat_pred, cross_covariance, cov_mat_cond)
  return(condmvn_output)
}

# Samples from the prior Gaussian Process
gp_prior_sample = function(N, x_pred, mean_func, kernel_func, ...) {
  # Samples N times from the prior GP at locations x_pred
  mean_vec = mean_func(x_pred, ...)
  cov_mat = kernel_func(x_pred, x_pred, ...)
  return(rmvn_svd(N, mean_vec, cov_mat))
}

# Samples from the posterior Gaussian Process
gp_posterior_sample = function(N, x_pred, x_cond, y_cond, mean_func, kernel_func, ...) {
  # Samples N from the posterior GP at locations x_pred
  # Conditioning f = y_cond at locations x_cond
  output = gp_posterior(x_pred, x_cond, y_cond, mean_func, kernel_func, ...)
  
  return(rmvn_svd(N, output$mean_vec, output$cov_mat))
}

# Evaluates the log-likelihood of the Gaussian process
gp_log_likelihood = function(x_cond, y_cond, mean_func, kernel_func, error_var = 0, ...) {
  # evaluate the log-likelihood of the Gaussian process
  mean_vector = mean_func(x_cond, ...)
  cov_mat = kernel_func(x_cond, x_cond, ...)
  
  n = nrow(cov_mat)
  error_cov_mat = diag(n) * error_var
  cov_mat = cov_mat + error_cov_mat
  
  log_likelihood = dmvn(y_cond, mean_vector, cov_mat)
  return(log_likelihood)
}