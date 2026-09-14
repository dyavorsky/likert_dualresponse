# gev_check.R -- does Lemma 1 survive when the inside-good errors are correlated?
#
# The paper assumes eta_itj iid Gumbel across alternatives. The natural
# relaxation is McFadden's GEV class, whose joint CDF is
#   F(e) = exp(-G(e^{-e_1}, ..., e^{-e_J})),  G homogeneous of degree 1.
#
# Claim (Appendix): under any GEV,
#   (i)  u* = max_j (V_j + e_j) ~ Gumbel(log G(e^V), 1);
#   (ii) argmax and max are independent.
#
# Verified here for a two-nest nested logit, which is GEV with
#   G(y) = sum_n ( sum_{j in nest n} y_j^{1/lambda_n} )^{lambda_n},
# and which is genuinely correlated whenever lambda_n < 1.
#
# Simulation uses the positive-stable mixture representation: conditional on
# S_n ~ PositiveStable(lambda_n), the errors in nest n are iid
# Gumbel(location = lambda_n * log S_n, scale = lambda_n).
#
# Run from repo root:  Rscript R/homogeneous/gev_check.R

set.seed(20260911)
N <- 2e6

# --- positive stable draws, Chambers-Mallows-Stuck ---------------------------
rposstable <- function(n, alpha) {
  U <- runif(n, 0, pi)
  W <- rexp(n)
  sin(alpha * U) / (sin(U))^(1 / alpha) *
    (sin((1 - alpha) * U) / W)^((1 - alpha) / alpha)
}

# --- design: 5 inside goods in 2 nests, deliberately unequal V ---------------
V      <- c(1.2, 0.4, -0.3, 0.9, -1.1)
nest   <- c(1, 1, 1, 2, 2)
lambda <- c(0.5, 0.8)          # dissimilarity; < 1 means correlated within nest

# --- GEV generator and its implied quantities --------------------------------
y <- exp(V)
Gval <- sum(sapply(1:2, function(n) sum(y[nest == n]^(1 / lambda[n]))^lambda[n]))

# GEV choice probabilities P_j = y_j G_j(y) / G(y), computed analytically
Gj <- numeric(length(V))
for (n in 1:2) {
  idx <- which(nest == n)
  inner <- sum(y[idx]^(1 / lambda[n]))
  Gj[idx] <- inner^(lambda[n] - 1) * y[idx]^(1 / lambda[n] - 1)
}
P_analytic <- y * Gj / Gval

# --- simulate ----------------------------------------------------------------
eta <- matrix(NA_real_, N, length(V))
for (n in 1:2) {
  idx <- which(nest == n)
  S   <- rposstable(N, lambda[n])
  loc <- lambda[n] * log(S)
  for (j in idx) {
    g <- -log(-log(runif(N)))                      # standard Gumbel
    eta[, j] <- loc + lambda[n] * g
  }
}
u     <- sweep(eta, 2, V, "+")
ustar <- apply(u, 1, max)
jstar <- max.col(u)

cat("=== (i) is max utility Gumbel(log G, 1)? ===\n")
cat(sprintf("  log G(e^V) analytic : %.4f\n", log(Gval)))
cat(sprintf("  location from mean  : %.4f   (mean - Euler)\n",
            mean(ustar) - 0.5772156649))
cat(sprintf("  sd of max           : %.4f   (Gumbel(.,1) sd = %.4f)\n",
            sd(ustar), pi / sqrt(6)))
ks <- suppressWarnings(ks.test(sample(ustar, 5000),
                               function(q) exp(-exp(-(q - log(Gval))))))
cat(sprintf("  KS vs Gumbel(logG,1): D = %.4f, p = %.3f\n",
            ks$statistic, ks$p.value))

cat("\n=== choice probabilities ===\n")
print(round(rbind(analytic = P_analytic,
                  simulated = as.numeric(table(jstar)) / N), 4))

cat("\n=== (ii) are argmax and max independent? ===\n")
cat("  mean of max, by chosen alternative (should all equal the overall mean):\n")
by_j <- tapply(ustar, jstar, mean)
cat(sprintf("   overall %.4f | %s\n", mean(ustar),
            paste(sprintf("j=%d: %.4f", seq_along(by_j), by_j), collapse = "  ")))
qs <- quantile(ustar, c(.1, .25, .5, .75, .9))
cat("  P(u* <= q | j*) vs P(u* <= q), at deciles/quartiles:\n")
for (q in qs) {
  cond <- tapply(ustar <= q, jstar, mean)
  cat(sprintf("   q=%6.3f  uncond %.4f | %s\n", q, mean(ustar <= q),
              paste(sprintf("%.4f", cond), collapse = " ")))
}
tab <- table(jstar, cut(ustar, breaks = c(-Inf, qs, Inf)))
cs  <- suppressWarnings(chisq.test(tab))
cat(sprintf("\n  chi-square test of independence: X2 = %.1f, df = %d, p = %.3f\n",
            cs$statistic, cs$parameter, cs$p.value))
cat("  (large p = cannot reject independence)\n")
