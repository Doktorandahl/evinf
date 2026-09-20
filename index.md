# evinf

`evinf` fits **extreme-value inflated** count-regression models: the
extreme-value and zero-inflated negative binomial (EVZINB) model and its
zero-inflation-free counterpart (EVINB), as described in

> Randahl, David, and Johan Vegelius. 2024. “Inference with Extremes:
> Accounting for Extreme Values in Count Regression Models.”
> *International Studies Quarterly* 68(4): sqae137.
> <https://doi.org/10.1093/isq/sqae137>

The package also provides bootstrap-based inference, prediction
(component probabilities, harmonic mean, quantiles), likelihood-ratio
tests for individual covariates, and comparison against plain
negative-binomial and ZINB models. It is designed to work with `broom` /
`generics` ([`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html)) and
`modelsummary`.

Since version 0.10.0 the EM tuning settings are bundled in
[`evinf_control()`](reference/evinf_control.md), the candidate range for
the extreme-value threshold `C_EV` is chosen from the data by default,
the count formula accepts an
[`offset()`](https://rdrr.io/r/stats/offset.html) term, and the fitted
models support the standard S3 methods
([`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`update()`](https://rdrr.io/r/stats/update.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)) as well as
[`predict_grid()`](reference/predict_grid.md),
[`marginal_effects()`](reference/marginal_effects.md) and
[`compare_fit()`](reference/compare_fit.md). See `NEWS.md` for the full
list.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("Doktorandahl/evinf")
```

## Usage

``` r

library(evinf)
data(genevzinb2)

future::plan(future::multisession, workers = 4) # bootstraps run in parallel

# c.lim pins the extreme-value threshold's candidate range (see c_profile());
# left NULL, it is chosen from the data.
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2,
                control = evinf_control(c.lim = c(173, 263)),
                n_bootstraps = 20)

summary(model)
tidy(model)                       # all components, one row per coefficient
glance(model)                     # goodness-of-fit
predict(model, type = "harmonic") # harmonic-mean prediction

lr_test(model, "x1")              # LR test for x1 across all components
comp <- compare_models(model)     # against NB and ZINB
compare_fit(comp)                 # paired bootstrap comparison of fit (AIC, BIC, ...)
```

Separate formulas can be given for the count, zero-inflation,
extreme-value inflation and Pareto components, e.g.

``` r

evzinb(y ~ x1 + x2, formula_zi = ~ x1, formula_pareto = ~ log1p(x3),
       data = genevzinb2)
```
