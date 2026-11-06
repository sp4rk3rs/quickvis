#### Load the package and run model checks when the suggested psych package is available
library(quickvis)
if (requireNamespace("psych", quietly = TRUE)) {

  #### Construct a known two-factor model with independently specified covariances
  lambda <- matrix(c(.8, .7, .6, 0, 0, 0, 0, 0, 0, .8, .7, .6), 6, 2)
  phi <- matrix(c(1, .3, .3, 1), 2, 2)
  sigma <- lambda %*% phi %*% t(lambda) + diag(1 - rowSums(lambda^2))
  dimnames(sigma) <- list(paste0("x", 1:6), paste0("x", 1:6))
  model <- "F1 =~ x1 + x2 + x3\nF2 =~ x4 + x5 + x6"
  set.seed(42)
  fit <- psych::CFA(model, sigma, n.obs = 500, n.iter = 2)

  #### Render and verify observed variance and covariance against the known model
  grDevices::png(tempfile(fileext = ".png"), width = 1600, height = 1200, res = 140)
  result <- path_diagram(fit, variance = "x1", covariance = c("x1", "x4"), tracing = TRUE)
  grDevices::dev.off()
  stopifnot(
    isTRUE(all.equal(unname(result$implied_covariance[paste0("x", 1:6), paste0("x", 1:6)]), unname(sigma))),
    abs(result$calculations$value[1] - 1) < 1e-8,
    abs(result$calculations$value[2] - .192) < 1e-8,
    abs(sum(result$traces$value[result$traces$calculation == "variance"]) - 1) < 1e-8,
    abs(sum(result$traces$value[result$traces$calculation == "covariance"]) - .192) < 1e-8
  )

  #### Check latent pairs, mixed pairs, single factors, and optional calculations
  grDevices::pdf(file = tempfile(fileext = ".pdf"), width = 12, height = 9)
  for (pair in list(c("F1", "F2"), c("F1", "x4"), c("x1", "F2"), c("x1", "x1"))) {
    value <- path_diagram(fit, covariance = pair, tracing = TRUE)
    stopifnot(abs(sum(value$traces$value) - value$calculations$value) < 1e-8)
  }
  plain <- path_diagram(fit)
  stopifnot(nrow(plain$calculations) == 0, nrow(plain$traces) == 0)
  latent_variance <- path_diagram(fit, variance = "F1", tracing = TRUE)
  stopifnot(abs(sum(latent_variance$traces$value) - 1) < 1e-8)
  comma_pair <- path_diagram(fit, covariance = "x1, x4")
  stopifnot(abs(comma_pair$calculations$value - .192) < 1e-8)
  single <- psych::CFA(sigma[1:3, 1:3], n.obs = 500, n.iter = 2)
  single_result <- path_diagram(single, variance = "x2", tracing = TRUE)
  stopifnot(abs(sum(single_result$traces$value) - 1) < 1e-8)

  #### Confirm negative single-factor loadings retain their signs
  signed_sigma <- sigma[1:3, 1:3] * outer(c(-1, 1, 1), c(-1, 1, 1))
  signed_fit <- psych::CFA(signed_sigma, n.obs = 500, n.iter = 2)
  signed_result <- path_diagram(signed_fit, covariance = c("x1", "x2"), tracing = TRUE)
  stopifnot(abs(signed_result$calculations$value + .56) < 1e-8,
            abs(sum(signed_result$traces$value) + .56) < 1e-8)

  #### Verify covariance-scale fitting and standardization with non-unit item variances
  scales <- c(2, 3, 4, 2, 3, 4)
  covariance_input <- sigma * outer(scales, scales)
  covariance_fit <- psych::CFA(model, covariance_input, n.obs = 500, n.iter = 2)
  unstandardized <- path_diagram(covariance_fit, std = "estimate", variance = "x1", covariance = c("x1", "x4"), tracing = TRUE)
  standardized <- path_diagram(covariance_fit, std = "all")
  latent_only <- path_diagram(covariance_fit, std = "lv", variance = "x1", tracing = TRUE)
  stopifnot(
    isTRUE(all.equal(unname(unstandardized$implied_covariance[paste0("x", 1:6), paste0("x", 1:6)]), unname(covariance_input))),
    isTRUE(all.equal(unname(standardized$implied_covariance[paste0("x", 1:6), paste0("x", 1:6)]), unname(sigma))),
    abs(unstandardized$calculations$value[1] - 4) < 1e-8,
    abs(unstandardized$calculations$value[2] - .768) < 1e-8,
    isTRUE(all.equal(latent_only$loadings, unstandardized$loadings)),
    isTRUE(all.equal(latent_only$residual_variance, unstandardized$residual_variance)),
    abs(latent_only$calculations$value - 4) < 1e-8,
    abs(sum(latent_only$traces$value) - 4) < 1e-8
  )

  #### Verify cross-loading variance terms include both symmetric contributions
  cross_fit <- fit
  cross_fit$cfa$Pattern[1, ] <- c(.6, .2)
  cross_fit$cfa$r <- cross_fit$cfa$Pattern %*% phi %*% t(cross_fit$cfa$Pattern) + diag(rep(.5, 6))
  cross_result <- path_diagram(cross_fit, variance = "x1", tracing = TRUE)
  stopifnot(abs(sum(cross_result$traces$value) - 1) < 1e-8)

  #### Check invalid arguments and restore the graphics device
  must_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
  must_error(path_diagram(fit, std = "estimate"))
  must_error(path_diagram(fit, tracing = TRUE))
  must_error(path_diagram(fit, digits = 1.5))
  must_error(path_diagram(fit, std = NA))
  must_error(path_diagram(fit, variance = "missing"))
  must_error(path_diagram(fit, covariance = c("x1", "missing")))
  must_error(path_diagram(list()))
  grDevices::dev.off()

  #### Report successful validation
  cat("All path diagram checks passed.\n")
}
