# quickvis

**quickvis** makes shaded loading tables for exploratory factor analysis (EFA)
and draws path diagrams for confirmatory factor analysis (CFA).

- `efa_table()` colors the loadings from a fitted EFA model.
- `compare_efa()` puts EFA solutions side by side on a shared color scale.
- `path_diagram()` plots a fitted `lavaan` CFA model with support for
  cross-loadings, bifactor layouts, and tracing of variances and covariances.

## Install

Clone or download this repository. From the directory containing the package
folder, run the following, using the folder name you downloaded if it differs
from `quickvis`.

```r
install.packages("remotes")
remotes::install_local("quickvis")
```

## Draw a CFA path diagram

```r
library(quickvis)

# Run once if needed: install.packages("lavaan")
model <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
'
fit <- lavaan::cfa(model, data = lavaan::HolzingerSwineford1939)
path_diagram(fit, covariance = c("x1", "x2"), tracing = TRUE)
```

The diagram supports single-group, single-level CFA models with continuous
indicators. See `?path_diagram` for model restrictions and plotting options.

## Compare EFA solutions

The [Harman example](examples/harman.R) uses `psych` to fit two-, three-, and
four-factor models and displays their loadings. The
[cross-loading and bifactor example](examples/path_diagram_cross_loadings.R)
and [Thurstone example](examples/thurstone_bifactor.R) show additional CFA
diagrams.

## License

[MIT](LICENSE.md)
