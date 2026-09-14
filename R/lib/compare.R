# compare.R -- model-comparison primitives shared by the aggregate and
# hierarchical estimators.

logmeanexp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

# Newton-Raftery: log p(y) = -logmeanexp(-ll_draws). trim = fraction of the
# smallest-likelihood draws dropped before averaging (stabilized variant).
lmd_nr <- function(ll_draws, trim = 0) {
  if (trim > 0) {
    keep <- ll_draws >= quantile(ll_draws, trim)
    ll_draws <- ll_draws[keep]
  }
  -logmeanexp(-ll_draws)
}
