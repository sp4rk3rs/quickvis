#### Load the package and run lavaan checks when the suggested dependency is available
library(quickvis)
if (requireNamespace("lavaan", quietly = TRUE)) {
  #### Reproduce the fixed factor variances and residual covariance in the reported model
  data <- lavaan::HolzingerSwineford1939
  observed <- c("visual", "cubes", "lozenge", "parag", "sentence", "word",
                "spd_add", "spd_cnt", "spd_dis")
  names(data)[match(paste0("x", 1:9), names(data))] <- observed
  model <- 'Visual_LF =~ NA*visual + cubes + lozenge
            Verbal =~ NA*parag + sentence + word
            Speed =~ NA*spd_add + spd_cnt + spd_dis
            Visual_LF ~~ 1*Visual_LF
            Verbal ~~ 1*Verbal
            Speed ~~ 1*Speed
            spd_add ~~ spd_cnt'
  fit <- lavaan::cfa(model, data = data, missing = "ML", estimator = "MLR")
  expected_covariance <- lavaan::lavInspect(fit, "cov.ov")

  #### Compare numerical matrices independently of lavaan's printing classes
  same_matrix <- function(a, b) {
    identical(dim(a), dim(b)) &&
      isTRUE(all.equal(as.numeric(a), as.numeric(b)))
  }

  #### Compare all three estimate scales with lavaan's independently standardized matrices
  grDevices::pdf(tempfile(fileext = ".pdf"), width = 13, height = 10)

  #### Check both documented covariance formats and explain misplaced variable names
  vector_pair <- path_diagram(fit, covariance = c("cubes", "word"), tracing = TRUE)
  string_pair <- path_diagram(fit, covariance = "cubes, word", tracing = TRUE)
  stopifnot(identical(vector_pair$calculations, string_pair$calculations))
  format_error <- tryCatch(path_diagram(fit, covariance = "cubes", "word", tracing = TRUE), error = identity)
  stopifnot(inherits(format_error, "error"), grepl("inside one argument", conditionMessage(format_error)))
  pair_error <- tryCatch(path_diagram(fit, covariance = "cubes"), error = identity)
  stopifnot(inherits(pair_error, "error"),
            grepl('covariance = c("cubes", "word")', conditionMessage(pair_error), fixed = TRUE))

  #### Check all numerical estimate scales after the syntax examples
  for (scale in c("all", "lv", "estimate")) {
    expected <- lavaan::lavInspect(fit, switch(scale, all = "std", lv = "std.lv", estimate = "est"))

    #### Use standardized residual covariances rather than residual correlations
    if (scale != "estimate") {
      solution <- lavaan::standardizedSolution(fit,
                                               type = paste0("std.", scale), cov.std = FALSE)
      rows <- which(solution$op == "~~" & solution$lhs %in% observed & solution$rhs %in% observed)
      for (row in rows) {
        a <- solution$lhs[row]
        b <- solution$rhs[row]
        expected$theta[a, b] <- expected$theta[b, a] <- solution$est.std[row]
      }
    }
    output <- path_diagram(fit, std = scale)
    target_covariance <- if (scale == "all") stats::cov2cor(expected_covariance) else expected_covariance
    stopifnot(
      identical(output$std, scale),
      same_matrix(output$loadings, expected$lambda),
      same_matrix(output$factor_covariance, expected$psi),
      same_matrix(output$residual_covariance, expected$theta),
      same_matrix(output$implied_covariance[observed, observed], target_covariance)
    )

    #### Verify all observed, latent, and mixed tracing totals against fitted covariances
    for (pair in list(c("spd_add", "spd_cnt"), c("spd_cnt", "spd_add"),
                      c("visual", "parag"), c("Visual_LF", "Verbal"),
                      c("Verbal", "spd_add"), c("spd_add", "spd_add"))) {
      traced <- path_diagram(fit, std = scale, covariance = pair, tracing = TRUE)
      stopifnot(abs(sum(traced$traces$value) - output$implied_covariance[pair[1], pair[2]]) < 1e-8)
    }
    variance <- path_diagram(fit, std = scale, variance = "spd_add", tracing = TRUE)
    stopifnot(abs(sum(variance$traces$value) - target_covariance["spd_add", "spd_add"]) < 1e-8)

    #### Ensure the correlated errors contribute their fitted covariance through a drawn edge
    residual_trace <- path_diagram(fit, std = scale,
                                   covariance = c("spd_add", "spd_cnt"), tracing = TRUE)
    residual_row <- grepl("<->", residual_trace$traces$path, fixed = TRUE)
    stopifnot(sum(residual_row) == 1L,
              abs(residual_trace$traces$value[residual_row] - expected$theta["spd_add", "spd_cnt"]) < 1e-8,
              any(grepl("^R", residual_trace$traces$edges[[which(residual_row)]])))
  }

  #### Verify factor standardization when marker identification leaves factor variances free
  marker_model <- 'Visual_LF =~ visual + cubes + lozenge
                   Verbal =~ parag + sentence + word
                   Speed =~ spd_add + spd_cnt + spd_dis
                   spd_add ~~ spd_cnt'
  marker_fit <- lavaan::cfa(marker_model, data = data, missing = "ML", estimator = "MLR")
  marker_expected <- lavaan::lavInspect(marker_fit, "std")
  marker_output <- path_diagram(marker_fit, variance = "Visual_LF", tracing = TRUE)
  marker_raw <- path_diagram(marker_fit, std = "estimate", variance = "Visual_LF", tracing = TRUE)
  marker_lv_expected <- lavaan::lavInspect(marker_fit, "std.lv")
  marker_lv <- path_diagram(marker_fit, std = "lv", variance = "Visual_LF", tracing = TRUE)
  stopifnot(same_matrix(marker_output$loadings, marker_expected$lambda),
            same_matrix(marker_output$factor_covariance, marker_expected$psi),
            abs(marker_output$calculations$value - 1) < 1e-8,
            abs(marker_raw$calculations$value - lavaan::lavInspect(marker_fit, "cov.lv")[1, 1]) < 1e-8,
            same_matrix(marker_lv$loadings, marker_lv_expected$lambda),
            same_matrix(marker_lv$factor_covariance, marker_lv_expected$psi),
            same_matrix(marker_lv$residual_covariance, marker_lv_expected$theta),
            same_matrix(marker_lv$implied_covariance[observed, observed],
                        lavaan::lavInspect(marker_fit, "cov.ov")),
            abs(marker_lv$calculations$value - 1) < 1e-8,
            abs(sum(marker_lv$traces$value) - 1) < 1e-8,
            !same_matrix(marker_lv$loadings, marker_raw$loadings),
            !same_matrix(marker_lv$loadings, marker_output$loadings))

  #### Confirm latent-only tracing for observed, mixed, and residual covariance requests
  for (pair in list(c("spd_cnt", "spd_cnt"), c("Visual_LF", "Verbal"),
                    c("Verbal", "spd_add"), c("spd_add", "spd_cnt"))) {
    traced <- path_diagram(marker_fit, std = "lv", covariance = pair, tracing = TRUE)
    stopifnot(abs(sum(traced$traces$value) - marker_lv$implied_covariance[pair[1], pair[2]]) < 1e-8)
  }

  #### Preserve legacy logical settings and reject ambiguous or invalid scales
  stopifnot(identical(path_diagram(marker_fit, std = TRUE), path_diagram(marker_fit, std = "all")),
            identical(path_diagram(marker_fit, std = FALSE), path_diagram(marker_fit, std = "estimate")))
  for (invalid in list("std.lv", "invalid", c("all", "lv"), character(), NA_character_, NA, 1, NULL)) {
    scale_error <- tryCatch(path_diagram(marker_fit, std = invalid), error = identity)
    stopifnot(inherits(scale_error, "error"), grepl('"all", "lv", or "estimate"', conditionMessage(scale_error)))
  }

  #### Confirm missing-data fitting does not substitute sample variances for fitted variances
  data$visual[1:10] <- NA_real_
  missing_fit <- lavaan::cfa(model, data = data, missing = "ML", estimator = "MLR")
  missing_output <- path_diagram(missing_fit, std = "estimate")
  stopifnot(same_matrix(missing_output$implied_covariance[observed, observed],
                       lavaan::lavInspect(missing_fit, "cov.ov")))

  #### Reject unsupported multi-group fits with an explicit explanation
  groups <- lavaan::cfa(model, data = data, group = "school", missing = "ML", estimator = "MLR")
  group_error <- tryCatch(path_diagram(groups), error = identity)
  stopifnot(inherits(group_error, "error"), grepl("single-group", conditionMessage(group_error)))
  grDevices::dev.off()

  #### Report successful validation
  cat("All lavaan path diagram checks passed.\n")
}
