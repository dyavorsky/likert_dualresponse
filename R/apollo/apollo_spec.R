# apollo_spec.R -- express the ordinal dual-response model in Apollo
# (Hess & Palma 2019) and fit it there.
#
# Purpose: an independent check of R/lib/likelihood.R. Apollo is the reference
# freeware toolkit in the choice-modelling community, so agreement between our
# estimates and Apollo's is stronger evidence than any internal test.
#
# The two models sit differently with respect to Apollo's built-in components:
#
#   Model B  ordered logit on the inclusive value. Native: apollo_mnl for the
#            forced choice, apollo_ol for the follow-up, joined by
#            apollo_combineModels. No custom code.
#
#   Model A  ordered complementary-log-log on the inclusive value. Apollo ships
#            apollo_ol (logistic) and apollo_op (normal) and nothing else, so
#            there is no native component. Supplied through apollo_ownModel,
#            which accepts a user-written likelihood and then hands it the full
#            Apollo machinery (estimation, s.e.'s, LR tests, prediction).
#
# The Model A probability below is deliberately written in its direct form,
#   Pr(y = w) = exp(-exp(-(c_w - mubar))) - exp(-exp(-(c_{w-1} - mubar))),
# rather than by reusing ord_prob_z() from R/lib/likelihood.R. Reusing our own
# code would make the comparison circular; the point is to check our
# tail-stable implementation against a transparent one.

if (!requireNamespace("apollo", quietly = TRUE)) {
  stop("The apollo package is required. Install with install.packages(\"apollo\").")
}
# Attached, not namespaced: Apollo parses the body of apollo_probabilities as
# source text to locate its model components, and a "apollo::" prefix defeats
# that parse (apollo_insertComponentName errors on the three-part call).
suppressPackageStartupMessages(library(apollo))

# ---- data ------------------------------------------------------------------

# Apollo stores one row per observation with every alternative's attributes on
# that row. Our design is long (n_tasks * J rows, task-major), so attribute p of
# alternative j for task t lives at row (t - 1) * J + j.
dual_to_apollo_db <- function(dat) {
  design <- dat$design
  n <- design$n_tasks
  J <- design$J
  X <- design$X
  attrs <- colnames(X)
  db <- data.frame(ID = seq_len(n), task = seq_len(n))
  for (j in seq_len(J)) {
    idx <- (seq_len(n) - 1L) * J + j
    for (p in attrs) db[[paste0(p, "_", j)]] <- X[idx, p]
  }
  db$choice <- dat$jstar
  db$y <- dat$y
  db
}

# ---- starting values -------------------------------------------------------

# Same rule as fit_dual_mle(): beta = 0, cut points backed out of the marginal
# category frequencies through the link. Kept identical so that any difference
# in the estimates is attributable to the likelihood or the optimiser, not to
# where each search started.
apollo_start <- function(dat) {
  design <- dat$design
  W <- dat$W
  attrs <- colnames(design$X)
  freq <- tabulate(dat$y, nbins = W)
  cumq <- pmin(pmax(cumsum(freq)[1:(W - 1)] / sum(freq), 1e-4), 1 - 1e-4)
  ginv <- if (dat$model == "B") qlogis else function(q) -log(-log(q))
  cut0 <- log(design$J) + ginv(cumq)
  if (W > 2) for (k in 2:(W - 1)) cut0[k] <- max(cut0[k], cut0[k - 1] + 1e-3)
  c(setNames(rep(0, length(attrs)), paste0("b_", attrs)),
    setNames(cut0, paste0("tau_", seq_len(W - 1))))
}

# ---- model definition ------------------------------------------------------

# Returns the apollo_probabilities function for one model. J, W and the
# attribute names are closed over, so the same builder serves any design.
make_apollo_probabilities <- function(dat) {
  J <- dat$design$J
  W <- dat$W
  attrs <- colnames(dat$design$X)
  model <- dat$model

  function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P <- list()

    # Forced choice: multinomial logit over the J displayed profiles. No
    # alternative-specific constants -- the cut points carry the only location.
    V <- list()
    for (j in seq_len(J)) {
      v <- 0
      for (p in attrs) v <- v + get(paste0("b_", p)) * get(paste0(p, "_", j))
      V[[paste0("alt", j)]] <- v
    }
    P[["choice"]] <- apollo_mnl(list(
      alternatives = setNames(seq_len(J), paste0("alt", seq_len(J))),
      avail        = 1,
      choiceVar    = choice,
      utilities    = V,
      componentName = "choice"
    ), functionality)

    # Inclusive value of the task. This is the single driver of the second
    # stage, and it carries no free loading parameter -- beta is shared with
    # the choice component. That restriction is the model's content.
    mubar <- log(Reduce(`+`, lapply(V, exp)))
    tau <- lapply(paste0("tau_", seq_len(W - 1)), function(z) get(z))

    if (model == "B") {
      P[["ord"]] <- apollo_ol(list(
        outcomeOrdered = y,
        utility        = mubar,
        tau            = tau,
        componentName  = "ord"
      ), functionality)
    } else {
      # Ordered cloglog, supplied by hand. Gumbel CDF of the standardised
      # bound, differenced across the two cut points bracketing the response.
      Fg <- function(z) exp(-exp(-z))
      lik <- function() {
        caug <- c(-Inf, unlist(tau, use.names = FALSE), Inf)
        hi <- caug[y + 1L] - mubar
        lo <- caug[y] - mubar
        Fg(hi) - Fg(lo)
      }
      P[["ord"]] <- apollo_ownModel(list(
        likelihood    = lik,
        componentName = "ord"
      ), functionality)
    }

    # Joint likelihood is the product of the two components (paper eq-joint).
    # No apollo_panelProd: the aggregate model treats tasks as independent, so
    # each task is its own ID.
    P <- apollo_combineModels(P, apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
}

# ---- fitting ---------------------------------------------------------------

# Fits one dataset in Apollo and returns the same shape fit_dual_mle() does.
fit_apollo <- function(dat, silent = TRUE) {
  database <<- dual_to_apollo_db(dat)
  apollo_initialise()

  apollo_control <<- list(
    modelName  = paste0("odr_model", dat$model),
    modelDescr = "Ordinal dual response",
    indivID    = "ID",
    outputDirectory = tempdir(),
    noValidation = TRUE,
    noDiagnostics = TRUE,
    panelData  = FALSE
  )
  apollo_beta <<- apollo_start(dat)
  apollo_fixed <<- c()
  apollo_probabilities <<- make_apollo_probabilities(dat)

  apollo_inputs <- apollo_validateInputs(silent = silent)
  est <- utils::capture.output(
    model <- apollo_estimate(
      apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
      estimate_settings = list(silent = silent, writeIter = FALSE,
                               hessianRoutine = "analytic")
    )
  )

  P <- length(colnames(dat$design$X))
  W <- dat$W
  bnm <- paste0("b_", colnames(dat$design$X))
  tnm <- paste0("tau_", seq_len(W - 1))
  list(
    beta = unname(model$estimate[bnm]),
    cut = unname(model$estimate[tnm]),
    se_beta = unname(model$se[bnm]),
    se_cut = unname(model$se[tnm]),
    ll = as.numeric(model$maximum),
    nll = -as.numeric(model$maximum),
    model = model
  )
}
