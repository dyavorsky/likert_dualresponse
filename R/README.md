# R/

Library code for the ordinal dual-response model. The **studies** that produce
the paper's results live in `../notebooks/`, not here: each notebook copies the
functions it needs from this directory, runs its experiment, and writes the
`R/output/*.rds` files the article reads.

## Layout

| directory | contents |
|---|---|
| `lib/` | shared machinery: `dgp.R` (designs, behavioral simulators), `likelihood.R` (joint likelihood, MLE fitter), `compare.R` (model-comparison primitives) |
| `homogeneous/` | aggregate-model pieces: `fisher.R` (exact information), `spec_tests.R` (the lambda and delta extensions), `bayes_agg.R` (aggregate posterior and marginal-likelihood estimators) |
| `hierarchical/` | panel pieces: `hb.R` (the sampler), `dgp_panel.R` (panel simulation), `compare_hb.R` (WAIC, holdout likelihood) |
| `apollo/` | the model expressed in Apollo, as an independent cross-check (Appendix E) |
| `output/` | results files, written by the notebooks and read by the article |

`homogeneous/gev_check.R` is the one remaining runnable script here. It verifies
the GEV lemma behind Appendix D numerically and is quoted in that appendix
rather than read from a results file.

## Where the studies went

Twelve runner scripts were retired once the notebooks took over as the source of
the article's numbers. Each notebook reproduces its script's results to
machine precision, so nothing was lost; the history is in git if a script is
ever wanted back.

| retired script | now in |
|---|---|
| `homogeneous/identification_test.R` | notebook 01 |
| `homogeneous/replication_study.R` | notebook 01 |
| `homogeneous/bias_scaling_check.R` | notebook 01 |
| `hierarchical/hb_coverage.R` | notebook 02 |
| `hierarchical/hb_validation.R` | notebook 02 |
| `homogeneous/dichotomization_study.R` | notebook 03 |
| `hierarchical/hb_simstudy.R` (cell `dichot`) | notebook 03 |
| `homogeneous/spec_test_study.R` | notebook 04 |
| `hierarchical/hb_simstudy.R` (five other cells) | notebook 05 |
| `homogeneous/design_guidance.R` | notebook 06 |
| `homogeneous/link_discrimination_study.R` | notebook 06 |
| `homogeneous/link_discrimination_squashed.R` | notebook 06 |
| `hierarchical/counterfactual_market.R` | notebook 07 |

## Rebuilding the results

`quarto render` from the repo root. The render order in `_quarto.yml` is
load-bearing: the notebooks run first and write `R/output/`, then `index.qmd`
reads it. Two further settings matter. `execute-dir: project` puts chunks in the
project root, and each notebook additionally anchors its output path on the
presence of `_quarto.yml`, because Quarto executes some passes from the file's
own directory regardless.

Editing a section under `sections/` does **not** invalidate the article's freeze,
since Quarto fingerprints `index.qmd` alone. After any such edit run
`rm -rf _freeze/index && quarto render` or the site will publish stale prose.
