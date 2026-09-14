# dgp.R -- design generation and behavioral simulation for the ordinal
# dual-response model. Simulation draws raw utilities and applies the two
# reporting rules directly (argmax for the choice; thresholded u* for Model A,
# thresholded u* - eta0 for Model B), NOT the derived likelihood -- so that
# parameter recovery via the closed-form likelihood also validates the
# derivation (factorization, link forms).

rgumbel <- function(n) -log(-log(runif(n)))

# Positive stable draws with index alpha in (0,1), Chambers-Mallows-Stuck.
# Used to simulate nested-logit (GEV) errors: conditional on S ~ PS(lambda),
# the errors within a nest are iid Gumbel(location = lambda*log S, scale =
# lambda), which reproduces the nested-logit joint distribution. Each marginal
# is still standard Gumbel, so nested and independent shocks are on the same
# scale and need no variance matching.
rposstable <- function(n, alpha) {
  U <- runif(n, 0, pi)
  W <- rexp(n)
  sin(alpha * U) / (sin(U))^(1 / alpha) *
    (sin((1 - alpha) * U) / W)^((1 - alpha) / alpha)
}

# Nested-logit shocks for an n x J matrix of tasks. `nest` is a list with
#   groups: length-J integer vector assigning each alternative to a nest
#   lambda: one dissimilarity per nest; lambda < 1 means correlated within.
rnested <- function(n, J, nest) {
  eta <- matrix(NA_real_, n, J)
  for (g in sort(unique(nest$groups))) {
    idx <- which(nest$groups == g)
    lam <- nest$lambda[g]
    loc <- lam * log(rposstable(n, lam))
    for (j in idx) eta[, j] <- loc + lam * rgumbel(n)
  }
  eta
}

# Design: two 3-level attributes (dummy-coded against omitted reference
# levels) plus one continuous price-like attribute. No intercept, no full
# dummy set -- satisfies the identification rank condition by construction.
# intercept = TRUE appends a column of ones (deliberate violation, for the
# negative control).
make_design <- function(n_tasks, J, seed, intercept = FALSE) {
  set.seed(seed)
  n_rows <- n_tasks * J
  a <- sample(1:3, n_rows, replace = TRUE)
  b <- sample(1:3, n_rows, replace = TRUE)
  price <- runif(n_rows, 0.5, 2.5)
  X <- cbind(
    a2 = as.numeric(a == 2),
    a3 = as.numeric(a == 3),
    b2 = as.numeric(b == 2),
    b3 = as.numeric(b == 3),
    price = price
  )
  if (intercept) X <- cbind(X, const = 1)
  list(X = X, n_tasks = n_tasks, J = J, P = ncol(X))
}

# The paper's identification condition (Proposition): [X, 1] must have full
# column rank P + 1.
id_rank_check <- function(X) {
  aug <- cbind(X, 1)
  rk <- qr(aug)$rank
  list(rank = rk, required = ncol(aug), pass = rk == ncol(aug))
}

# ---------------------------------------------------------------------------
# Panel (respondent x task) structure for hierarchical models
# ---------------------------------------------------------------------------

# Panel design: N respondents x T tasks x J alternatives, same attribute
# scheme as make_design. Rows are respondent-major, task-major, alt-fastest.
make_panel_design <- function(N, T_tasks, J, seed, intercept = FALSE) {
  des <- make_design(N * T_tasks, J, seed = seed, intercept = intercept)
  des$N <- N
  des$T_tasks <- T_tasks
  des$resp <- rep(seq_len(N), each = T_tasks)          # respondent of each task
  des$resp_row <- rep(des$resp, each = J)              # respondent of each X row
  des
}

# Simulate the behavioral process. `cut` is on the utility scale for both
# models (for Model A, cut = -log(-log(alpha)) i.e. alpha = exp(-exp(-cut))).
#   model2: "inclusive" (the paper's second stage) or "2max" (DR-2Max-style:
#           second-stage latent is V_{j*} + fresh logistic noise)
#   lambda: slope on the inclusive value in the second-stage latent (Model B
#           semantics; lambda = 1 is the theory's value)
# The model2 and lambda branches mirror simulate_dual_hb() exactly, including
# the order in which random numbers are drawn, so the two simulators agree
# draw for draw on a single-respondent panel.
simulate_dual <- function(design, beta, cut, model = c("A", "B"), seed,
                          model2 = c("inclusive", "2max"), lambda = 1) {
  model <- match.arg(model)
  model2 <- match.arg(model2)
  stopifnot(length(beta) == design$P, !is.unsorted(cut, strictly = TRUE))
  set.seed(seed)
  n <- design$n_tasks
  J <- design$J
  V <- matrix(design$X %*% beta, nrow = n, ncol = J, byrow = TRUE)
  u <- V + matrix(rgumbel(n * J), n, J)
  jstar <- max.col(u, ties.method = "first")
  ustar <- u[cbind(seq_len(n), jstar)]

  m <- do.call(pmax, as.data.frame(V))
  logS <- m + log(rowSums(exp(V - m)))

  latent <- if (model2 == "2max") {
    V[cbind(seq_len(n), jstar)] + rlogis(n)
  } else if (lambda != 1) {
    rlogis(n, location = lambda * logS)
  } else if (model == "A") {
    ustar
  } else {
    ustar - rgumbel(n)
  }
  y <- findInterval(latent, cut) + 1L
  list(
    design = design, jstar = jstar, y = y, W = length(cut) + 1L,
    model = model, beta_true = beta, cut_true = cut,
    model2 = model2, lambda = lambda
  )
}
