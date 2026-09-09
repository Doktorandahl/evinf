# evinf

<!-- badges: start -->
[![R-CMD-check](https://github.com/Doktorandahl/evinf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Doktorandahl/evinf/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/Doktorandahl/evinf/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/Doktorandahl/evinf/actions/workflows/test-coverage.yaml)
[![Codecov test coverage](https://codecov.io/gh/Doktorandahl/evinf/graph/badge.svg)](https://app.codecov.io/gh/Doktorandahl/evinf)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

`evinf` fits **extreme-value inflated** count-regression models: the
extreme-value and zero-inflated negative binomial (EVZINB) model and its
zero-inflation-free counterpart (EVINB), as described in

> Randahl, David, and Johan Vegelius. 2024. "Inference with Extremes:
> Accounting for Extreme Values in Count Regression Models."
> *International Studies Quarterly* 68(4): sqae137.
> <https://doi.org/10.1093/isq/sqae137>

The package also provides bootstrap-based inference, prediction (component
probabilities, harmonic mean, quantiles), likelihood-ratio tests for individual
covariates, and comparison against plain negative-binomial and ZINB models. It
is designed to work with `broom` / `generics` (`tidy()`, `glance()`) and
`modelsummary`.

Since version 0.10.0 the EM tuning settings are bundled in `evinf_control()`,
the candidate range for the extreme-value threshold `C_EV` is chosen from the
data by default, the count formula accepts an `offset()` term, and the fitted
models support the standard S3 methods (`coef()`, `vcov()`, `confint()`,
`AIC()`, `predict()`, `residuals()`, `simulate()`, `update()`, `plot()`) as
well as `predict_grid()`, `marginal_effects()` and `compare_fit()`. See
`NEWS.md` for the full list.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Doktorandahl/evinf")
```

## Usage

```r
library(evinf)
data(genevzinb2)

# Fit an EVZINB model with 100 bootstrap replicates
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 100)

summary(model)
tidy(model)                       # all components, one row per coefficient
glance(model)                     # goodness-of-fit
predict(model, type = "harmonic") # harmonic-mean prediction

lr_test(model, "x1")              # LR test for x1 across all components
compare_models(model)             # against NB and ZINB
```

Separate formulas can be given for the count, zero-inflation, extreme-value
inflation and Pareto components, e.g.

```r
evzinb(y ~ x1 + x2, formula_zi = ~ x1, formula_pareto = ~ log1p(x3),
       data = genevzinb2)
```
