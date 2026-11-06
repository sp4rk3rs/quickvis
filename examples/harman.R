library(quickvis)

#### Check the packages needed to run this example
if (!requireNamespace("psych", quietly = TRUE)) {
  stop(
    "Install the psych package before running this example.",
    call. = FALSE
  )
}

#### Load the example correlation matrix included with R
data("Harman74.cor", package = "datasets")

#### Fit two-, three-, and four-factor solutions
harman_models <- lapply(2:4, function(n) {
  psych::fa(
    r = Harman74.cor$cov,
    n.obs = Harman74.cor$n.obs,
    nfactors = n,
    fm = "minres",
    rotate = "oblimin"
  )
})

#### Display the solutions with shared shading and space between models
print(efa_table(harman_models[[1]]))
print(
  compare_efa(
    "Two factors" = harman_models[[1]],
    "Three factors" = harman_models[[2]],
    "Four factors" = harman_models[[3]],
    model_gap = 24
  )
)
