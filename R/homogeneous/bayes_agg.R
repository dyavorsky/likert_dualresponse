# bayes_agg.R -- Bayesian estimation and marginal-likelihood estimators for
# the AGGREGATE model. Requires R/lib/likelihood.R and R/lib/compare.R.


# --- aggregate-model posterior + bridge sampling -----------------------------

# RW-MH sampler for the aggregate model's posterior. Prior: N(0, tau2 * I) on
# the working parameters (beta, c1, log-diffs). Proposal covariance is the
# MLE's inverse Hessian, scaled by 2.4^2/d (adaptive during burn-in).
fit_dual_bayes_agg <- function(dat, iter = 20000, burn = 5000, thin = 5,
                               tau2 = 100, seed = NULL, verbose = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  mle <- fit_dual_mle(dat)
  d <- length(mle$par)
  H <- optimHess(mle$par, function(p) negloglik(p, dat))
  Vprop <- chol2inv(chol(H))
  Lp <- t(chol(Vprop)) * (2.4 / sqrt(d))
  lpost <- function(p) -negloglik(p, dat) - 0.5 * sum(p^2) / tau2
  cur <- mle$par
  cur_lp <- lpost(cur)
  nkeep <- floor((iter - burn) / thin)
  draws <- matrix(NA_real_, nkeep, d)
  lp_draws <- numeric(nkeep)
  acc <- 0; ki <- 0; sc <- 1
  for (it in seq_len(iter)) {
    prop <- cur + as.numeric(Lp %*% rnorm(d)) * sc
    prop_lp <- lpost(prop)
    if (log(runif(1)) < prop_lp - cur_lp) {
      cur <- prop; cur_lp <- prop_lp; acc <- acc + 1
    }
    if (it <= burn && it %% 200 == 0) {
      sc <- min(max(sc * exp(1.5 * (acc / 200 - 0.234)), 0.1), 10)
      acc <- 0
    }
    if (it == burn) acc <- 0
    if (it > burn && (it - burn) %% thin == 0) {
      ki <- ki + 1
      draws[ki, ] <- cur
      lp_draws[ki] <- cur_lp
    }
  }
  list(draws = draws, lp = lp_draws, accept = acc / (iter - burn),
       mle = mle, tau2 = tau2, dat = dat)
}

# Bridge-sampling LMD for the aggregate posterior (needs `bridgesampling`).
bridge_lmd_agg <- function(fit_agg) {
  if (!requireNamespace("bridgesampling", quietly = TRUE)) {
    warning("bridgesampling not installed; returning NA")
    return(NA_real_)
  }
  dat <- fit_agg$dat
  tau2 <- fit_agg$tau2
  lpost_fn <- function(pars, data) {
    -negloglik(pars, data) - 0.5 * sum(pars^2) / tau2 -
      0.5 * length(pars) * log(2 * pi * tau2)   # include prior normalizing constant
  }
  samples <- fit_agg$draws
  colnames(samples) <- paste0("p", seq_len(ncol(samples)))
  lb <- rep(-Inf, ncol(samples)); ub <- rep(Inf, ncol(samples))
  names(lb) <- names(ub) <- colnames(samples)
  bs <- bridgesampling::bridge_sampler(
    samples = samples, log_posterior = lpost_fn, data = dat,
    lb = lb, ub = ub, silent = TRUE
  )
  bs$logml
}

# NR estimate on the aggregate posterior draws (for the NR-vs-bridge check).
lmd_nr_agg <- function(fit_agg, trim = 0) {
  ll <- apply(fit_agg$draws, 1, function(p) -negloglik(p, fit_agg$dat))
  lmd_nr(ll, trim = trim)
}

# Laplace approximation to the aggregate LMD (independent cross-check on the
# bridge estimate): ll(theta_hat) + log prior(theta_hat) + (d/2) log 2pi
# - 0.5 log |H|, H = observed information at the MLE (posterior mode under a
# flat-enough prior).
lmd_laplace_agg <- function(fit_agg) {
  mle <- fit_agg$mle
  d <- length(mle$par)
  H <- optimHess(mle$par, function(p) negloglik(p, fit_agg$dat))
  lprior <- -0.5 * sum(mle$par^2) / fit_agg$tau2 - 0.5 * d * log(2 * pi * fit_agg$tau2)
  -mle$nll + lprior + 0.5 * d * log(2 * pi) - 0.5 * determinant(H, logarithm = TRUE)$modulus
}
