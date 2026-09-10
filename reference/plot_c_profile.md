# Plot the log-likelihood profile over the candidate values of C_EV

Plot the log-likelihood profile over the candidate values of C_EV

## Usage

``` r
plot_c_profile(object)
```

## Arguments

- object:

  A fitted `evzinb` / `evinb` object.

## Value

A `ggplot` object.

## Examples

``` r
# \donttest{
data(genevzinb2)
model <- evzinb(y ~ x1 + x2 + x3, data = genevzinb2, n_bootstraps = 5)
#> evinf: using a data-driven candidate range for C_EV: [173, 263]. Pass `c.lim` / `control = evinf_control(c.lim = ...)` to override.
plot_c_profile(model)

# }
```
