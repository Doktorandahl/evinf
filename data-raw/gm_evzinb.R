## Regenerates data/gm_evzinb.rda
##
## gm_evzinb is the `gof_map` passed to modelsummary so that the goodness-of-fit
## rows produced by glance.evzinb() / glance.evinb() are labelled and formatted
## sensibly. It must stay in sync with the columns those glance methods return.

gm_evzinb <- tibble::tribble(
  ~raw,                  ~clean,                 ~fmt,
  "parameter",           "C_EV",                 0,
  "alpha",               "alpha_nb",             2,
  "nobs",                "obs",                  0,
  "npar",                "par",                  0,
  "converged",           "converged",            0,
  "n_bootstraps",        "n_bootstraps",         0,
  "n_failed_bootstraps", "n_failed_bootstraps",  0,
  "aic",                 "AIC",                  1,
  "bic",                 "BIC",                  1,
  "logLik",              "logLik",               2
)

usethis::use_data(gm_evzinb, overwrite = TRUE)
