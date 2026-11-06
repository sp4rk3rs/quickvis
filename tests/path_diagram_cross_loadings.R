#### Test fitted cross-loading and bifactor diagrams against known population models
library(quickvis)
if (requireNamespace("lavaan", quietly = TRUE)) {
  same_matrix <- function(a, b) {
    identical(dim(a), dim(b)) && isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-6))
  }
  must_error <- function(expr, message) {
    error <- tryCatch(force(expr), error = identity)
    stopifnot(inherits(error, "error"), grepl(message, conditionMessage(error), fixed = TRUE))
  }
  fit_population <- function(model, lambda, phi, theta) {
    sigma <- lambda %*% phi %*% t(lambda) + theta
    dimnames(sigma) <- list(rownames(lambda), rownames(lambda))
    lavaan::cfa(model, sample.cov = sigma, sample.nobs = 10000,
                sample.cov.rescale = FALSE, std.lv = TRUE)
  }
  loading_edges <- function(output) Filter(function(edge) edge$type == "loading", output$edges)
  check_tracing <- function(output) {
    stopifnot(all(unlist(output$traces$edges) %in% names(output$edges)))
    for (kind in output$calculations$calculation) {
      stopifnot(abs(sum(output$traces$value[output$traces$calculation == kind]) -
                      output$calculations$value[output$calculations$calculation == kind]) < 1e-8)
    }
  }
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 12, height = 9)

  #### A single observed variable loads on two correlated factors
  lambda <- cbind(F1 = c(.8, .7, .6, 0, 0, 0), F2 = c(0, 0, .3, .8, .7, .6))
  rownames(lambda) <- paste0("x", 1:6)
  phi <- matrix(c(1, .3, .3, 1), 2)
  theta <- diag(rep(.5, 6))
  cross <- fit_population("F1 =~ x1 + x2 + x3\nF2 =~ x3 + x4 + x5 + x6", lambda, phi, theta)
  for (scale in c("all", "lv", "estimate")) {
    output <- path_diagram(cross, std = scale, variance = "x3",
                           covariance = c("x3", "x4"), tracing = TRUE)
    expected <- lavaan::lavInspect(cross, switch(scale, all = "std", lv = "std.lv", estimate = "est"))
    stopifnot(output$layout == "tree", same_matrix(output$loadings, expected$lambda),
              sum(output$nodes$name == "x3") == 1L, length(loading_edges(output)) == 7L,
              identical(vapply(output$edges[c("L3:1", "L3:2")], `[[`, "", "to"),
                        c("L3:1" = "x3", "L3:2" = "x3")),
              output$edges[["L3:1"]]$x[2] != output$edges[["L3:2"]]$x[2])
    check_tracing(output)
    cross_term <- output$traces$calculation == "variance" & grepl("2 * ", output$traces$expression, fixed = TRUE)
    stopifnot(sum(cross_term) == 1L,
              abs(output$traces$value[cross_term] -
                    2 * output$loadings["x3", 1] * output$loadings["x3", 2] *
                    output$factor_covariance[1, 2]) < 1e-8)
    if (scale == "estimate") stopifnot(abs(output$traces$value[cross_term] - .108) < 1e-6)
    reverse <- path_diagram(cross, std = scale, covariance = c("x4", "x3"), tracing = TRUE)
    check_tracing(reverse)
    stopifnot(abs(reverse$calculations$value - output$calculations$value[2]) < 1e-8)
  }

  #### A negative cross-loading preserves its sign and the signed cross term
  lambda["x3", "F2"] <- -.3
  negative <- fit_population("F1 =~ x1 + x2 + x3\nF2 =~ .8*x4 + x5 + x6 + x3", lambda, phi, theta)
  output <- path_diagram(negative, std = "estimate", variance = "x3", tracing = TRUE)
  check_tracing(output)
  stopifnot(output$edges[["L3:2"]]$value < 0,
            any(output$traces$value < 0))

  #### A general factor appears below indicators even when its loadings are strongest
  lambda <- cbind(g = rep(.7, 9), S1 = c(.5, .4, .6, rep(0, 6)),
                  S2 = c(rep(0, 3), .5, .4, .6, rep(0, 3)),
                  S3 = c(rep(0, 6), .5, .4, .6))
  rownames(lambda) <- paste0("x", 1:9)
  model <- 'g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
            S1 =~ x1 + x2 + x3
            S2 =~ x4 + x5 + x6
            S3 =~ x7 + x8 + x9
            g ~~ 0*S1 + 0*S2 + 0*S3'
  phi <- diag(4)
  phi[2, 3] <- phi[3, 2] <- .2
  phi[2, 4] <- phi[4, 2] <- .1
  phi[3, 4] <- phi[4, 3] <- .15
  theta <- diag(rep(.5, 9))
  theta[1, 2] <- theta[2, 1] <- .05
  bifactor <- fit_population(paste(model, "x1 ~~ x2", sep = "\n"), lambda, phi, theta)
  for (scale in c("all", "lv", "estimate")) {
    output <- path_diagram(bifactor, std = scale, variance = "x1",
                           covariance = c("x1", "x2"), tracing = TRUE)
    nodes <- output$nodes
    g <- nodes[nodes$name == "g", ]
    items <- nodes[nodes$type == "observed", ]
    stopifnot(output$layout == "bifactor", nrow(nodes) == 13L,
              length(loading_edges(output)) == 18L,
              g$y < min(items$y), all(diff(items$x) > 0),
              all(nodes$y[nodes$name %in% c("S1", "S2", "S3")] > max(items$y)),
              all(vapply(loading_edges(output)[paste0("L", 1:9, ":1")],
                         function(edge) diff(edge$y) > 0, logical(1))),
              all(c("C2:3", "C2:4", "C3:4", "R1:2") %in% names(output$edges)))
    check_tracing(output)
    expected <- lavaan::lavInspect(bifactor, "cov.ov")
    if (scale == "all") expected <- stats::cov2cor(expected)
    stopifnot(same_matrix(output$implied_covariance[rownames(lambda), rownames(lambda)], expected))
    mixed <- path_diagram(bifactor, std = scale, covariance = c("g", "x1"), tracing = TRUE)
    check_tracing(mixed)
    tree <- path_diagram(bifactor, std = scale, layout = "tree")
    stopifnot(tree$layout == "tree", length(loading_edges(tree)) == 18L,
              all(tree$nodes$y[tree$nodes$type == "latent"] > max(tree$nodes$y[tree$nodes$type == "observed"])),
              identical(tree$implied_covariance, output$implied_covariance))
  }

  #### Retain nonzero general/specific covariances and allow explicit placement
  phi[1, 2:4] <- phi[2:4, 1] <- c(.1, -.1, .05)
  correlated_model <- sub("g ~~ 0*S1 + 0*S2 + 0*S3", "g ~~ .1*S1 + -.1*S2 + .05*S3", model, fixed = TRUE)
  correlated <- fit_population(correlated_model, lambda, phi, diag(rep(.5, 9)))
  output <- path_diagram(correlated, general = "g", variance = "x1", tracing = TRUE)
  check_tracing(output)
  stopifnot(all(c("C1:2", "C1:3", "C1:4") %in% names(output$edges)))

  #### Match the sketch with two upper factors and one lower general factor
  two_lambda <- lambda[1:6, 1:3]
  two_phi <- diag(3)
  two_phi[2, 3] <- two_phi[3, 2] <- .2
  two_model <- 'g =~ .7*x1 + .7*x2 + .7*x3 + .7*x4 + .7*x5 + .7*x6
                S1 =~ x1 + x2 + x3
                S2 =~ x4 + x5 + x6
                g ~~ 0*S1 + 0*S2'
  two_bifactor <- fit_population(two_model, two_lambda, two_phi, diag(rep(.5, 6)))
  output <- path_diagram(two_bifactor, covariance = c("x1", "x4"), tracing = TRUE)
  stopifnot(output$layout == "bifactor", nrow(output$nodes) == 9L,
            length(loading_edges(output)) == 12L, "C2:3" %in% names(output$edges))
  check_tracing(output)

  #### Ambiguous general factors require a name; partial general factors can be named
  ambiguous_model <- 'A =~ .8*x1 + .7*x2 + .6*x3
                       B =~ .5*x1 + .4*x2 + .3*x3
                       A ~~ 0*B'
  lambda <- cbind(A = c(.8, .7, .6), B = c(.5, .4, .3))
  rownames(lambda) <- paste0("x", 1:3)
  ambiguous <- fit_population(ambiguous_model, lambda, diag(2), diag(rep(.5, 3)))
  stopifnot(path_diagram(ambiguous)$layout == "tree",
            path_diagram(ambiguous, layout = "bifactor", general = "B")$layout == "bifactor")
  must_error(path_diagram(ambiguous, layout = "bifactor"), "unique general factor")
  stopifnot(path_diagram(cross, general = "F1")$layout == "bifactor")
  must_error(path_diagram(cross, general = "missing"), "Available factors")
  must_error(path_diagram(cross, general = c("F1", "F2")), "one latent factor")
  must_error(path_diagram(cross, general = "F1", layout = "tree"), "Use `layout")

  #### Distinguish a freely estimated zero cross-loading from a fixed absent path
  zero_model <- 'F1 =~ x1 + x2 + x3
                 F2 =~ cross*x3 + x4 + x5 + x6
                 cross == 0'
  lambda <- cbind(F1 = c(.8, .7, .6, 0, 0, 0), F2 = c(0, 0, 0, .8, .7, .6))
  rownames(lambda) <- paste0("x", 1:6)
  zero <- fit_population(zero_model, lambda, diag(2), diag(rep(.5, 6)))
  output <- path_diagram(zero)
  stopifnot("L3:2" %in% names(output$edges), abs(output$edges[["L3:2"]]$value) < 1e-8)
  fixed_zero <- fit_population(sub("cross*x3", "0*x3", sub("cross == 0", "", zero_model, fixed = TRUE), fixed = TRUE),
                               lambda, diag(2), diag(rep(.5, 6)))
  stopifnot(!"L3:2" %in% names(path_diagram(fixed_zero)$edges))
  grDevices::dev.off()
  cat("All cross-loading and bifactor path diagram checks passed.\n")
}
