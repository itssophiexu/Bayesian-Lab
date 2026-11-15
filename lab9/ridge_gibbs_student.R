# Complete exactly the five right-hand-side expressions marked below.
# The model uses an independent fixed-variance normal prior for the slopes.
# The intercept is included in the coefficient block and is not penalized.

ridge_gibbs <- function(X, y, v, n_iter = 3000L, burn = 1000L,
                        a0 = 2, b0 = 1, initial_s2 = 1, seed = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  validate_gibbs_inputs(X, y, v, n_iter, burn, a0, b0, initial_s2)
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  Z <- cbind("(Intercept)" = 1, X)
  q <- ncol(Z)
  G <- crossprod(Z)
  h <- drop(crossprod(Z, y))

  P <- diag(c(0, rep(1 / v, p)))
  require_completed(P, "prior precision")
  stopifnot(is.matrix(P), all(dim(P) == c(q, q)))

  theta <- c(mean(y), rep(0, p))
  s2 <- initial_s2
  theta_trace <- matrix(NA_real_, n_iter + 1L, q)
  sigma2_trace <- numeric(n_iter + 1L)
  theta_trace[1, ] <- theta
  sigma2_trace[1] <- s2
  colnames(theta_trace) <- colnames(Z)
  for (iteration in seq_len(n_iter)) {
    A <- G / s2 + P
    h_cond <- h / s2
    require_completed(A, "conditional precision")
    require_completed(h_cond, "conditional mean numerator")
    stopifnot(is.matrix(A), all(dim(A) == c(q, q)), length(h_cond) == q)
    theta <- draw_from_precision(A, h_cond)
    rss <- sum((y - drop(Z %*% theta))^2)
    shape <- a0 + n / 2
    rate <- b0 + rss / 2
    require_completed(shape, "variance shape")
    require_completed(rate, "variance rate")
    stopifnot(length(shape) == 1L, shape > 0, length(rate) == 1L, rate > 0)
    s2 <- 1 / rgamma(1, shape = shape, rate = rate)
    theta_trace[iteration + 1L, ] <- theta
    sigma2_trace[iteration + 1L] <- s2
  }
  retained <- seq.int(burn + 2L, n_iter + 1L)
  list(theta = theta_trace[retained, , drop = FALSE],
       sigma2 = sigma2_trace[retained], theta_trace = theta_trace,
       sigma2_trace = sigma2_trace,
       settings = list(v = v, a0 = a0, b0 = b0, n_iter = n_iter,
                       burn = burn, thinning = 1L, initial_s2 = initial_s2,
                       seed = seed))
}
