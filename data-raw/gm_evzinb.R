## Regenerates data/gm_evzinb.rda
##
## gm_evzinb is the `gof_map` passed to modelsummary so that the goodness-of-fit
## rows produced by glance.evzinb() / glance.evinb() are labelled and formatted
## sensibly. It is exactly gof_map_evinf(); see R/gof_map_evinf.R. Keep this in
## sync with the columns those glance methods return.

devtools::load_all(".")

gm_evzinb <- gof_map_evinf()

usethis::use_data(gm_evzinb, overwrite = TRUE)
