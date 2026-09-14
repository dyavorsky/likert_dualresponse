# compare_apollo.R -- cross-validation of R/lib/likelihood.R against Apollo.
#
# Fits the same simulated datasets two ways: with our maximum-likelihood
# routine, and with the same model expressed in Apollo (Hess & Palma 2019).
# Agreement is evidence that the likelihood in the paper is implemented
# correctly, checked against the reference toolkit of the choice-modelling
# community rather than against ourselves.
#
# Three configurations:
#   B / W = 5   ordered logit on the inclusive value -- native apollo_ol
#   A / W = 5   ordered cloglog                      -- custom apollo_ownModel
#   B / W = 2   the binary dual-response special case (the incumbent model)
#
# Run from repo root: Rscript R/apollo/compare_apollo.R
# Saves R/output/apollo_compare.rds

source(file.path("R", "lib", "dgp.R"))
source(file.path("R", "lib", "likelihood.R"))
source(file.path("R", "apollo", "apollo_spec.R"))
dir.create(file.path("R", "output"), showWarnings = FALSE)

R_REPS <- 5
N_TASKS <- 2000
J <- 4

beta_true <- c(a2 = 0.8, a3 = -0.5, b2 = 0.4, b3 = 1.0, price = -0.9)
alpha_to_cut <- function(a) -log(-log(a))

configs <- list(
  list(key = "B_W5", model = "B", cut = c(-1.0, 0.2, 1.2, 2.2)),
  list(key = "A_W5", model = "A", cut = alpha_to_cut(c(0.10, 0.30, 0.60, 0.85))),
  list(key = "B_W2", model = "B", cut = 0.5)
)

run_config <- function(cfg) {
  d <- lapply(seq_len(R_REPS), function(r) {
    des <- make_design(N_TASKS, J, seed = 4000 + r)
    dat <- simulate_dual(des, beta_true, cfg$cut, model = cfg$model,
                         seed = 5000 + r)
    fo <- fit_dual_mle(dat)
    fa <- fit_apollo(dat)
    c(beta = max(abs(fo$beta - fa$beta)),
      cut = max(abs(fo$cut - fa$cut)),
      se_beta = max(abs(fo$se_beta - fa$se_beta)),
      se_cut = max(abs(fo$se_cut - fa$se_cut)),
      ll = abs(-fo$nll - fa$ll))
  })
  m <- do.call(rbind, d)
  data.frame(config = cfg$key, model = cfg$model, W = length(cfg$cut) + 1L,
             component = if (cfg$model == "A") "apollo_ownModel" else "apollo_ol",
             max_d_beta = max(m[, "beta"]), max_d_cut = max(m[, "cut"]),
             max_d_se_beta = max(m[, "se_beta"]), max_d_se_cut = max(m[, "se_cut"]),
             max_d_ll = max(m[, "ll"]))
}

cat(sprintf("Apollo cross-validation: %d replications per configuration, n = %d tasks\n",
            R_REPS, N_TASKS))
cat(sprintf("apollo version %s\n\n", as.character(packageVersion("apollo"))))

res <- do.call(rbind, lapply(configs, run_config))
rownames(res) <- NULL

print(res, digits = 3)

worst <- max(res$max_d_beta, res$max_d_cut)
cat(sprintf("\nLargest disagreement on any point estimate: %.2e\n", worst))
cat(sprintf("Largest disagreement on any standard error:  %.2e\n",
            max(res$max_d_se_beta, res$max_d_se_cut)))
cat(sprintf("Largest disagreement on the log-likelihood:  %.2e\n", max(res$max_d_ll)))

out <- list(res = res, settings = list(R_REPS = R_REPS, N_TASKS = N_TASKS, J = J,
                                       apollo_version = as.character(packageVersion("apollo"))))
saveRDS(out, file.path("R", "output", "apollo_compare.rds"))
cat("Saved: R/output/apollo_compare.rds\n")
