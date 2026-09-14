# likelihood.R -- joint log-likelihood of the ordinal dual-response model and
# an MLE fitter with identification diagnostics.
#
# Joint per-task likelihood (paper eq-joint):
#   [ exp(V_j*) / S ] * [ G(c_y - mubar) - G(c_{y-1} - mubar) ]
# with S = sum_j exp(V_j), mubar = log S, and link G:
#   Model A: standard Gumbel CDF  G(z) = exp(-exp(-z))   (probability report)
#   Model B: standard logistic    G(z) = plogis(z)       (graded comparison)
#
# Parameterization for unconstrained optimization:
#   par = ( beta[1..P], c_1, log(c_2 - c_1), ..., log(c_{W-1} - c_{W-2}) )

par_to_cut <- function(par, P, W) {
  if (W == 2) return(par[P + 1])
  par[P + 1] + c(0, cumsum(exp(par[(P + 2):(P + W - 1)])))
}

cut_to_par <- function(cut) {
  if (length(cut) == 1) return(cut)
  c(cut[1], log(diff(cut)))
}

row_max <- function(M) do.call(pmax, as.data.frame(M))

# Interval probability G(hi) - G(lo) with tail-stable formulas, from the
# standardized bounds directly (lo = c_{w-1} - mubar, hi = c_w - mubar).
ord_prob_z <- function(lo, hi, model) {
  if (model == "B") {
    # plogis(hi) - plogis(lo) = plogis(hi) * plogis(-lo) * (1 - exp(lo - hi))
    plogis(hi) * plogis(-lo) * (-expm1(lo - hi))
  } else {
    # exp(-e^{-hi}) - exp(-e^{-lo}) = exp(-e^{-hi}) * (1 - exp(-(e^{-lo} - e^{-hi})))
    ea <- exp(-lo)
    eb <- exp(-hi)
    p <- exp(-eb) * (-expm1(-(ea - eb)))
    # deep left tail: exp(-hi) overflows and both CDFs underflow to zero
    p[!is.finite(eb)] <- 0
    p
  }
}

# P( y = w | mubar ) with common cut points; c_0 = -Inf, c_W = +Inf.
ord_prob <- function(mubar, cut, y, model) {
  caug <- c(-Inf, cut, Inf)
  ord_prob_z(caug[y] - mubar, caug[y + 1L] - mubar, model)
}

negloglik <- function(par, dat) {
  design <- dat$design
  P <- design$P
  W <- dat$W
  n <- design$n_tasks
  J <- design$J
  beta <- par[1:P]
  cut <- par_to_cut(par, P, W)

  V <- matrix(design$X %*% beta, nrow = n, ncol = J, byrow = TRUE)
  m <- row_max(V)
  logS <- m + log(rowSums(exp(V - m)))

  ll_choice <- V[cbind(seq_len(n), dat$jstar)] - logS
  p_ord <- ord_prob(logS, cut, dat$y, dat$model)
  ll_ord <- log(pmax(p_ord, 1e-312))

  -(sum(ll_choice) + sum(ll_ord))
}

num_grad <- function(f, x, eps = 1e-6) {
  vapply(seq_along(x), function(k) {
    h <- eps * max(1, abs(x[k]))
    xp <- x; xp[k] <- xp[k] + h
    xm <- x; xm[k] <- xm[k] - h
    (f(xp) - f(xm)) / (2 * h)
  }, numeric(1))
}

# MLE with observed-information diagnostics. Starting values: beta = 0 (so
# mubar = log J for every task) and cut points backed out from the marginal
# category frequencies through the link.
fit_dual_mle <- function(dat, start = NULL) {
  design <- dat$design
  P <- design$P
  W <- dat$W

  if (is.null(start)) {
    freq <- tabulate(dat$y, nbins = W)
    cumq <- pmin(pmax(cumsum(freq)[1:(W - 1)] / sum(freq), 1e-4), 1 - 1e-4)
    ginv <- if (dat$model == "B") qlogis else function(q) -log(-log(q))
    cut0 <- log(design$J) + ginv(cumq)
    if (W > 2) {
      for (k in 2:(W - 1)) cut0[k] <- max(cut0[k], cut0[k - 1] + 1e-3)
    }
    start <- c(rep(0, P), cut_to_par(cut0))
  }

  fn <- function(p) negloglik(p, dat)
  opt <- optim(start, fn, method = "BFGS",
               control = list(maxit = 1000, reltol = 1e-12))
  H <- optimHess(opt$par, fn)
  g <- num_grad(fn, opt$par)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values

  cut_hat <- par_to_cut(opt$par, P, W)
  # Delta method for the cut points: dc_w/dc_1 = 1; dc_w/d theta_m = d_m for m <= w.
  Jc <- matrix(0, W - 1, length(opt$par))
  Jc[, P + 1] <- 1
  if (W > 2) {
    d <- exp(opt$par[(P + 2):(P + W - 1)])
    for (w in 2:(W - 1)) Jc[w, (P + 2):(P + w)] <- d[1:(w - 1)]
  }
  Vpar <- tryCatch(solve(H), error = function(e) matrix(NA, length(opt$par), length(opt$par)))
  se_beta <- sqrt(diag(Vpar)[1:P])
  se_cut <- sqrt(diag(Jc %*% Vpar %*% t(Jc)))

  out <- list(
    beta = opt$par[1:P], se_beta = se_beta,
    cut = cut_hat, se_cut = se_cut,
    nll = opt$value, convergence = opt$convergence,
    grad_max = max(abs(g)), eigen = ev,
    min_eig = min(ev), cond = max(ev) / max(min(ev), .Machine$double.eps),
    par = opt$par, counts = opt$counts
  )
  if (dat$model == "A") {
    alpha_hat <- exp(-exp(-cut_hat))
    # d alpha / d c = alpha * (-log alpha)
    out$alpha <- alpha_hat
    out$se_alpha <- se_cut * alpha_hat * (-log(alpha_hat))
  }
  out
}
