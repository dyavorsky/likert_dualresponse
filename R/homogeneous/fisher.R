# fisher.R -- exact Fisher information for the ordinal dual-response model,
# implementing (and generalizing to the full (beta, cut) block) the paper's
# information decomposition: per task,
#   I = [ MNL block ] + [ ordinal-stage block ],
# where the beta-beta part of the ordinal block is omega(mubar; c) * xbar xbar'
# with omega = sum_w (g_w - g_{w-1})^2 / pi_w  (paper eq-fisher).
#
# Parameter order: (beta[1..P], c_1, ..., c_{W-1}).

link_cdf <- function(z, model) if (model == "B") plogis(z) else exp(-exp(-z))
link_pdf <- function(z, model) {
  if (model == "B") dlogis(z) else exp(-(z + exp(-z)))
}

# Information contribution of one task with design rows Xt (J x P).
fisher_task <- function(Xt, beta, cut, model) {
  P <- ncol(Xt)
  W <- length(cut) + 1L
  K <- P + W - 1L
  V <- as.numeric(Xt %*% beta)
  eV <- exp(V - max(V))
  p <- eV / sum(eV)
  mubar <- max(V) + log(sum(eV))

  I <- matrix(0, K, K)
  # MNL block (choice among inside goods)
  I[1:P, 1:P] <- t(Xt) %*% (diag(p) - tcrossprod(p)) %*% Xt

  # Ordinal-stage block over the full (beta, cut) vector
  z <- cut - mubar
  Gz <- c(0, link_cdf(z, model), 1)
  gz <- c(0, link_pdf(z, model), 0)
  pi_w <- diff(Gz)
  xbar <- as.numeric(t(p) %*% Xt)
  for (w in seq_len(W)) {
    s <- numeric(K)
    s[1:P] <- -(gz[w + 1] - gz[w]) / pi_w[w] * xbar
    if (w <= W - 1) s[P + w] <- s[P + w] + gz[w + 1] / pi_w[w]
    if (w >= 2)     s[P + w - 1] <- s[P + w - 1] - gz[w] / pi_w[w]
    I <- I + pi_w[w] * tcrossprod(s)
  }
  I
}

# Total information over a design (rows task-major, J per task).
fisher_total <- function(design, beta, cut, model) {
  n <- design$n_tasks
  J <- design$J
  K <- design$P + length(cut)
  I <- matrix(0, K, K)
  for (t in seq_len(n)) {
    rows <- ((t - 1) * J + 1):(t * J)
    I <- I + fisher_task(design$X[rows, , drop = FALSE], beta, cut, model)
  }
  I
}

# Second-stage-only information weight omega(mubar; cut) per task (the scalar
# in eq-fisher). Useful for quantifying how much of the second response's
# information dichotomization discards, independent of the design.
omega_weight <- function(mubar, cut, model) {
  z <- outer(mubar, cut, function(m, c) c - m)
  Gz <- cbind(0, link_cdf(z, model), 1)
  gz <- cbind(0, link_pdf(z, model), 0)
  W <- length(cut) + 1L
  om <- 0
  for (w in seq_len(W)) {
    om <- om + (gz[, w + 1] - gz[, w])^2 / (Gz[, w + 1] - Gz[, w])
  }
  om
}

# Asymptotic variance of the predicted exceedance probability
# h = P(y >= k | task) = 1 - G(c_{k-1} - mubar) for a single evaluation task,
# given a parameter covariance Vcov whose cut block corresponds to `cut`
# (positions cut_pos, a vector of the indices of c_1.. within the parameter).
pred_var_exceed <- function(Xt, beta, cut, k, model, Vcov, P) {
  V <- as.numeric(Xt %*% beta)
  m <- max(V)
  mubar <- m + log(sum(exp(V - m)))
  p <- exp(V - mubar)
  xbar <- as.numeric(t(p) %*% Xt)
  g <- link_pdf(cut[k - 1] - mubar, model)
  grad <- numeric(nrow(Vcov))
  grad[1:P] <- g * xbar
  grad[P + k - 1] <- -g
  as.numeric(t(grad) %*% Vcov %*% grad)
}

# True exceedance probability for a task.
true_exceed <- function(Xt, beta, cut, k, model) {
  V <- as.numeric(Xt %*% beta)
  m <- max(V)
  mubar <- m + log(sum(exp(V - m)))
  1 - link_cdf(cut[k - 1] - mubar, model)
}
