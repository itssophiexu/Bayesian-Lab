# STA 402, Week of November 16, 2026
# Start R in the lab folder, then source("lab_starter.R").
# Sections 1 and 2 run immediately. Complete the five TODOs in sections 3
# and 4. Then change FALSE to TRUE in section 5 and source this file again.

# 1. Read and inspect the observed counts.
counts <- read.csv("data/nfl_2026_field_goals.csv")
band_score <- c("30-39" = -1, "40-49" = 0, "50-59" = 1)
stopifnot(
  all(counts$distance_band %in% names(band_score)),
  all(counts$attempts >= counts$made & counts$made >= 0),
  !anyDuplicated(counts[c("kicker", "distance_band")])
)
counts <- counts[counts$attempts > 0, ]
kicker_names <- sort(unique(counts$kicker))
focus <- c("Jake Bates", "Spencer Shrader")
band_totals <- aggregate(cbind(made, attempts) ~ distance_band, counts, sum)
band_totals$raw_rate <- with(band_totals, made / attempts)
long_counts <- counts[counts$distance_band == "50-59", c("kicker", "made", "attempts")]
long_counts$raw_rate <- with(long_counts, made / attempts)
print(band_totals, row.names = FALSE)
print(long_counts[match(focus, long_counts$kicker), ], row.names = FALSE)

# 2. Expand each kicker-by-band cell to binary observations.
# The within-cell row order is arbitrary and carries no time information.
cell <- rep(seq_len(nrow(counts)), counts$attempts)
kicks <- data.frame(
  kicker = counts$kicker[cell],
  group = match(counts$kicker[cell], kicker_names),
  x = unname(band_score[counts$distance_band[cell]]),
  y = as.integer(sequence(counts$attempts) <= counts$made[cell])
)
stopifnot(nrow(kicks) == sum(counts$attempts), sum(kicks$y) == sum(counts$made))

# 3. Gibbs sampler. Normal distributions in the handout use variances,
# while R's rnorm uses standard deviations.
draw_latent <- function(eta, y) {
  sign_y <- 2 * y - 1
  signed_eta <- sign_y * eta
  # sign_y * z is Normal(signed_eta, 1) restricted to positive values.
  # After centering, the standard Normal upper tail at -signed_eta
  # has probability pnorm(signed_eta). Sample uniformly on that tail
  # probability scale. Work in logs for numerical stability.
  log_tail <- log(runif(length(y))) + pnorm(signed_eta, log.p = TRUE)

  # TODO 1: Invert log_tail using qnorm with lower.tail = FALSE and
  # log.p = TRUE. Restore the mean, then restore z's sign.
  z <- stop("Complete TODO 1: truncated Normal draw")
  stopifnot(all(is.finite(z)), all(sign_y * z > 0))
  z
}

gibbs_probit <- function(b, seed, mu_start, beta_start, tau2_start,
                         iterations = 4000L, burn = 1000L) {
  stopifnot(iterations > burn, burn >= 0, b > 0, tau2_start > 0)
  set.seed(seed)
  g <- kicks$group
  x <- kicks$x
  y <- kicks$y
  J <- length(kicker_names)
  n_j <- tabulate(g, nbins = J)
  mu <- mu_start
  beta <- beta_start
  tau2 <- tau2_start
  alpha <- rep(mu, J)
  draws <- matrix(NA_real_, nrow = iterations - burn, ncol = J + 3L,
                  dimnames = list(NULL, c("mu", "beta", "tau2", kicker_names)))

  for (iteration in seq_len(iterations)) {
    z <- draw_latent(alpha[g] + beta * x, y)

    # TODO 2: Complete the vector of conditional variances and means for alpha.
    # rowsum adds the working responses within kicker in increasing group order.
    sum_by_kicker <- as.numeric(rowsum(z - beta * x, g, reorder = TRUE))
    alpha_var <- stop("Complete TODO 2: kicker variances")
    alpha_mean <- stop("Complete TODO 2: kicker means")
    alpha <- rnorm(J, mean = alpha_mean, sd = sqrt(alpha_var))

    # TODO 3: Complete the variance and mean for the shared distance effect.
    beta_var <- stop("Complete TODO 3: distance effect variance")
    beta_mean <- stop("Complete TODO 3: distance effect mean")
    beta <- rnorm(1, mean = beta_mean, sd = sqrt(beta_var))

    # Given alpha and tau2, combine J Normal terms with the Normal(1, 1) prior.
    mu_var <- 1 / (J / tau2 + 1)
    mu_mean <- mu_var * (sum(alpha) / tau2 + 1)
    mu <- rnorm(1, mean = mu_mean, sd = sqrt(mu_var))

    # TODO 4: Complete the shape and scale of the inverse-gamma update.
    # Under the handout's convention, 1 / tau2 has a Gamma conditional
    # with rate equal to the inverse-gamma scale.
    tau_shape <- stop("Complete TODO 4: inverse-gamma shape")
    tau_scale <- stop("Complete TODO 4: inverse-gamma scale")
    tau2 <- 1 / rgamma(1, shape = tau_shape, rate = tau_scale)

    if (iteration > burn) {
      draws[iteration - burn, ] <- c(mu, beta, tau2, alpha)
    }
  }
  draws
}

fit_prior <- function(b, seed) {
  list(
    gibbs_probit(b, seed, mu_start = 0.3, beta_start = 0.5, tau2_start = 0.03),
    gibbs_probit(b, seed + 1L, mu_start = 1.7, beta_start = -0.5, tau2_start = 0.5)
  )
}

plot_traces <- function(chains) {
  saved_par <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  on.exit(par(saved_par))
  trace_names <- c("mu", "beta", "tau2", focus[1])
  trace_labels <- c("mu", "beta", "tau squared", paste("alpha:", focus[1]))
  for (j in seq_along(trace_names)) {
    name <- trace_names[j]
    trace <- cbind(chains[[1]][, name], chains[[2]][, name])
    matplot(trace, type = "l", lty = 1, col = c("#143D59", "#B45C31"),
            xlab = "Retained draw", ylab = trace_labels[j])
    if (j == 1L) {
      legend("topright", c("Chain 1", "Chain 2"), lty = 1,
             col = c("#143D59", "#B45C31"), bty = "n")
    }
  }
}

# 4. Summarize 50-59-yard probability draws and predict a 10-attempt count.
long_probabilities <- function(chains, kicker) {
  draws <- rbind(chains[[1]], chains[[2]])
  pnorm(draws[, kicker] + draws[, "beta"])
}

summarize_long <- function(chains, prior, kicker) {
  p <- long_probabilities(chains, kicker)
  interval <- quantile(p, c(0.05, 0.95), names = FALSE)
  observed <- long_counts[match(kicker, long_counts$kicker), ]
  data.frame(prior = prior, kicker = kicker, made = observed$made,
             attempts = observed$attempts, raw_rate = observed$raw_rate,
             posterior_mean = mean(p), lower_90 = interval[1],
             upper_90 = interval[2])
}

predict_ten <- function(chains, prior, kicker) {
  p <- long_probabilities(chains, kicker)
  # TODO 5: For every draw of p, draw the count of makes in 10 new attempts.
  makes <- stop("Complete TODO 5: posterior predictive count")
  interval <- quantile(makes, c(0.05, 0.95), names = FALSE, type = 1)
  data.frame(prior = prior, kicker = kicker, expected_makes = mean(makes),
             lower_90 = interval[1], upper_90 = interval[2],
             prob_at_least_7 = mean(makes >= 7))
}

# 5. Run the analysis after filling the five TODOs.
if (FALSE) {
  fits <- list(strong = fit_prior(b = 0.05, seed = 40201L),
               weak = fit_prior(b = 0.50, seed = 40203L))

  # Look at both chains within each fit before interpreting these summaries.
  # Run plot_traces(fits$strong) and plot_traces(fits$weak) in the console.
  parameter_summary <- do.call(rbind, lapply(names(fits), function(prior) {
    draws <- rbind(fits[[prior]][[1]], fits[[prior]][[2]])
    data.frame(prior = prior, median_tau2 = median(draws[, "tau2"]),
               prob_beta_negative = mean(draws[, "beta"] < 0))
  }))
  probability_summary <- do.call(rbind, lapply(names(fits), function(prior) {
    do.call(rbind, lapply(focus, function(kicker) {
      summarize_long(fits[[prior]], prior, kicker)
    }))
  }))
  set.seed(40205L)
  predictive_summary <- do.call(rbind, lapply(names(fits), function(prior) {
    do.call(rbind, lapply(focus, function(kicker) {
      predict_ten(fits[[prior]], prior, kicker)
    }))
  }))
  print(parameter_summary, digits = 3, row.names = FALSE)
  print(probability_summary, digits = 3, row.names = FALSE)
  print(predictive_summary, digits = 3, row.names = FALSE)
}
