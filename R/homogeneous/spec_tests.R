# spec_tests.R -- the paper's two sharp specification tests, implemented as
# one-parameter extensions of the aggregate likelihood:
#
#   lambda (unit slope): P(y <= w) = G(c_w - lambda * mubar).
#     Theory fixes lambda = 1 (both responses driven by the same utility
#     draw); a free slope is identified by how the ordinal distribution
#     shifts as task attractiveness varies. Test: H0 lambda = 1.
#
#   delta (conditioning): second-stage predictor mubar + delta*(V_{j*} - mubar).
#     delta = 0 is our inclusive-value model; delta = 1 is the DR-2Max
#     conditioning (second stage reads the chosen option's utility).
#     Tests: H0 delta = 0 (our model) and H0 delta = 1 (DR-2Max).
#
# Both are Wald + LR tests from the aggregate MLE, with the extended fits
# warm-started at the restricted solution.
#
# Requires: dgp.R, likelihood.R.

# Extended negative log-likelihood. ext = "lambda" or "delta"; if `fix` is
# NULL the extension parameter is the last element of `par`, otherwise it is
# held at `fix` and `par` has the standard layout (beta, zeta).
negloglik_ext <- function(par, dat, ext = c("lambda", "delta"), fix = NULL) {
  ext <- match.arg(ext)
  design <- dat$design
  P <- design$P
  W <- dat$W
  n <- design$n_tasks
  J <- design$J
  if (is.null(fix)) {
    theta <- par[P + W]
    par_core <- par[1:(P + W - 1)]
  } else {
    theta <- fix
    par_core <- par
  }
  beta <- par_core[1:P]
  cut <- par_to_cut(par_core, P, W)

  V <- matrix(design$X %*% beta, nrow = n, ncol = J, byrow = TRUE)
  m <- row_max(V)
  logS <- m + log(rowSums(exp(V - m)))
  Vstar <- V[cbind(seq_len(n), dat$jstar)]
  ll_choice <- Vstar - logS

  pred <- if (ext == "lambda") theta * logS else logS + theta * (Vstar - logS)
  p_ord <- ord_prob(pred, cut, dat$y, dat$model)
  -(sum(ll_choice) + sum(log(pmax(p_ord, 1e-312))))
}

# Fit the extended model, warm-started at the restricted MLE.
fit_dual_mle_ext <- function(dat, ext, start_core, theta0) {
  fn <- function(p) negloglik_ext(p, dat, ext = ext)
  opt <- optim(c(start_core, theta0), fn, method = "BFGS",
               control = list(maxit = 1000, reltol = 1e-12))
  H <- optimHess(opt$par, fn)
  Vc <- tryCatch(solve(H), error = function(e) matrix(NA, length(opt$par), length(opt$par)))
  k <- length(opt$par)
  list(par = opt$par, theta = opt$par[k], se_theta = sqrt(Vc[k, k]),
       nll = opt$value, convergence = opt$convergence)
}

# Run the full test battery on one aggregate dataset. Returns a data frame of
# tests: estimate, se, Wald z, LR stat, p-values.
spec_test_battery <- function(dat, mle = NULL) {
  if (is.null(mle)) mle <- fit_dual_mle(dat)
  start_core <- mle$par

  # lambda test (H0: lambda = 1)
  fl <- fit_dual_mle_ext(dat, "lambda", start_core, theta0 = 1)
  lr_l <- 2 * (mle$nll - fl$nll)
  # delta test (H0: delta = 0); restricted-at-0 fit IS the standard MLE
  fd <- fit_dual_mle_ext(dat, "delta", start_core, theta0 = 0)
  lr_d0 <- 2 * (mle$nll - fd$nll)
  # H0: delta = 1 (DR-2Max) -- constrained fit with delta fixed at 1
  fn_d1 <- function(p) negloglik_ext(p, dat, ext = "delta", fix = 1)
  opt_d1 <- optim(start_core, fn_d1, method = "BFGS",
                  control = list(maxit = 1000, reltol = 1e-12))
  lr_d1 <- 2 * (opt_d1$value - fd$nll)

  data.frame(
    test = c("lambda = 1", "delta = 0", "delta = 1 (DR-2Max)"),
    estimate = c(fl$theta, fd$theta, fd$theta),
    se = c(fl$se_theta, fd$se_theta, fd$se_theta),
    wald_z = c((fl$theta - 1) / fl$se_theta,
               fd$theta / fd$se_theta,
               (fd$theta - 1) / fd$se_theta),
    lr = c(lr_l, lr_d0, lr_d1),
    p_wald = 2 * pnorm(-abs(c((fl$theta - 1) / fl$se_theta,
                              fd$theta / fd$se_theta,
                              (fd$theta - 1) / fd$se_theta))),
    p_lr = pchisq(c(lr_l, lr_d0, lr_d1), df = 1, lower.tail = FALSE)
  )
}
