# compare_hb.R -- model comparison for the hierarchical model: WAIC from
# per-respondent pointwise log-likelihood draws, and the log posterior-
# predictive density of held-out tasks. Requires R/lib/compare.R and
# R/hierarchical/hb.R.


# WAIC from an (ndraws x N) matrix of per-respondent log-likelihoods.
waic <- function(ll_resp_draws) {
  lppd <- sum(apply(ll_resp_draws, 2, logmeanexp))
  p_waic <- sum(apply(ll_resp_draws, 2, var))
  list(waic = -2 * (lppd - p_waic), lppd = lppd, p_waic = p_waic,
       elpd = lppd - p_waic)
}

# Log posterior-predictive density of holdout tasks, summed over respondents:
# sum_i log mean_r exp( ll_i(holdout | phi_i^(r), cut^(r)) ).
holdout_ll <- function(fit, dat_holdout) {
  stopifnot(!is.null(fit$draws$Phi))
  P <- fit$settings$P
  W <- fit$settings$W
  het_cut <- fit$settings$het_cut
  nkeep <- dim(fit$draws$Phi)[1]
  N <- dim(fit$draws$Phi)[2]
  ll <- matrix(NA_real_, nkeep, N)
  for (r in seq_len(nkeep)) {
    Phi_r <- fit$draws$Phi[r, , , drop = TRUE]
    if (is.null(dim(Phi_r))) Phi_r <- matrix(Phi_r, N, 1)
    cutmat <- if (het_cut) phi_to_cutmat(Phi_r, P, W)
              else matrix(fit$draws$cut[r, ], N, W - 1, byrow = TRUE)
    ll[r, ] <- resp_loglik_all(dat_holdout, Phi_r[, 1:P, drop = FALSE], cutmat,
                               fit$settings$choice_only)
  }
  list(total = sum(apply(ll, 2, logmeanexp)),
       by_resp = apply(ll, 2, logmeanexp))
}
