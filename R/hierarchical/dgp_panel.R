# dgp_panel.R -- panel (hierarchical) simulation: individual coefficient
# draws, the behavioral process on a panel design, and the estimation /
# holdout split. Requires R/lib/dgp.R for designs and shock primitives.


# Draw individual coefficients from the population model.
draw_betas <- function(N, beta_bar, Sigma, seed) {
  set.seed(seed)
  P <- length(beta_bar)
  Z <- matrix(rnorm(N * P), N, P)
  sweep(Z %*% chol(Sigma), 2, beta_bar, "+")
}

# Simulate the behavioral process on a panel design.
#   Bmat:   N x P matrix of individual coefficients
#   cut:    common cut vector (length W-1) or N x (W-1) matrix (heterogeneous)
#   model:  "A" (probability report), "B" (graded comparison)
#   model2: "inclusive" (the paper's second stage) or "2max" (DR-2Max-style:
#           second-stage latent is V_{j*} + fresh logistic noise)
#   lambda: slope on the inclusive value in the second-stage latent (Model B
#           semantics; lambda = 1 is the theory's value)
#   shock:  "gumbel" or "normal" (variance-matched, for misspecification runs)
#   nest:   NULL, or list(groups = <length-J nest ids>, lambda = <one per
#           nest>) to draw correlated nested-logit (GEV) first-stage
#           shocks instead of independent ones
#   drift:  NULL, list(type = "linear", delta = d) shifting effective cut
#           points by d * (t - 1) across a respondent's task sequence, or
#           list(type = "update", rho = r) for Hung-style probabilistic
#           updating of the outside-good location to the previous task's
#           inclusive value with probability plogis(rho)
simulate_dual_hb <- function(design, Bmat, cut, model = c("A", "B"), seed,
                             model2 = c("inclusive", "2max"), lambda = 1,
                             shock = c("gumbel", "normal"), drift = NULL,
                             nest = NULL) {
  model <- match.arg(model)
  model2 <- match.arg(model2)
  shock <- match.arg(shock)
  set.seed(seed)
  n <- design$n_tasks
  J <- design$J
  N <- design$N
  stopifnot(nrow(Bmat) == N, ncol(Bmat) == design$P)

  cutmat <- if (is.matrix(cut)) cut else matrix(cut, N, length(cut), byrow = TRUE)
  W <- ncol(cutmat) + 1L

  Bx <- Bmat[design$resp_row, , drop = FALSE]
  V <- matrix(rowSums(design$X * Bx), nrow = n, ncol = J, byrow = TRUE)
  u <- if (!is.null(nest)) {
    V + rnested(n, J, nest)
  } else {
    eta <- if (shock == "gumbel") rgumbel(n * J) else rnorm(n * J, sd = pi / sqrt(6))
    V + matrix(eta, n, J)
  }
  jstar <- max.col(u, ties.method = "first")
  ustar <- u[cbind(seq_len(n), jstar)]

  m <- do.call(pmax, as.data.frame(V))
  logS <- m + log(rowSums(exp(V - m)))

  latent <- if (model2 == "2max") {
    V[cbind(seq_len(n), jstar)] + rlogis(n)
  } else if (lambda != 1) {
    rlogis(n, location = lambda * logS)      # statistical alternative for the slope test
  } else if (model == "A") {
    ustar
  } else {
    ustar - rgumbel(n)
  }

  shift <- numeric(n)
  if (!is.null(drift)) {
    T_tasks <- design$T_tasks
    tt <- rep(seq_len(T_tasks), times = N)
    if (drift$type == "linear") {
      shift <- drift$delta * (tt - 1)
    } else if (drift$type == "update") {
      psi <- numeric(n)
      for (i in seq_len(N)) {
        rows <- which(design$resp == i)
        for (k in seq_along(rows)[-1]) {
          psi[rows[k]] <- if (runif(1) < plogis(drift$rho)) logS[rows[k - 1]] else psi[rows[k - 1]]
        }
      }
      shift <- psi
    }
  }

  eff_cut <- cutmat[design$resp, , drop = FALSE] + shift
  y <- 1L + rowSums(latent >= eff_cut)

  list(design = design, resp = design$resp, resp_row = design$resp_row,
       N = N, jstar = jstar, y = as.integer(y), W = W, model = model,
       Bmat_true = Bmat, cut_true = cut, model2 = model2, lambda = lambda,
       shock = shock, drift = drift, nest = nest)
}

# Split a panel dataset into estimation and holdout tasks per respondent
# (the last `n_holdout` tasks of each respondent are held out).
split_holdout <- function(dat, n_holdout) {
  des <- dat$design
  T_tasks <- des$T_tasks
  keep_t <- rep(seq_len(T_tasks) <= T_tasks - n_holdout, times = des$N)
  subset_panel <- function(keep) {
    rows <- rep(keep, each = des$J)
    d <- des
    d$X <- des$X[rows, , drop = FALSE]
    d$n_tasks <- sum(keep)
    d$T_tasks <- d$n_tasks / des$N
    d$resp <- des$resp[keep]
    d$resp_row <- rep(d$resp, each = des$J)
    out <- dat
    out$design <- d
    out$resp <- d$resp
    out$resp_row <- d$resp_row
    out$jstar <- dat$jstar[keep]
    out$y <- dat$y[keep]
    out
  }
  list(est = subset_panel(keep_t), holdout = subset_panel(!keep_t))
}
