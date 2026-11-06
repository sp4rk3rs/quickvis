#' Create a CFA path diagram with cross-loadings or a general factor
#'
#' @param object A fitted model returned by [lavaan::cfa()].
#' @param digits Number of decimal places to display.
#' @param std Estimate scale: `"all"` (default) standardizes latent and observed
#'   variables, `"lv"` standardizes latent variables only, and `"estimate"` uses
#'   unstandardized estimates. These correspond to lavaan's `Std.all`, `Std.lv`,
#'   and `Estimate` columns. For compatibility, `TRUE` means `"all"` and
#'   `FALSE` means `"estimate"`.
#' @param variance One quoted latent or observed variable name, for example
#'   `variance = "Visual_LF"`. Defaults to `NULL` (no variance request).
#' @param covariance Two quoted variable names inside a character vector,
#'   for example `covariance = c("cubes", "word")`. Alternatively, use one
#'   comma-separated string: `covariance = "cubes, word"`. Both names belong
#'   inside this one argument; `covariance = "cubes", "word"` is not valid.
#'   Defaults to `NULL` (no covariance request).
#' @param tracing Logical. Whether to highlight the paths contributing
#'   to the requested variance or covariance using colored arrows
#'   and an explanatory legend.
#' @param layout Diagram arrangement: `"auto"` (default) places a general factor
#'   below the indicators when a bifactor pattern is detected; otherwise factors
#'   appear above the indicators. Use `"tree"` to keep all factors above, or
#'   `"bifactor"` to request the general-factor arrangement explicitly.
#' @param general Optional name of the general factor to place below the
#'   indicators. With `NULL`, a unique factor loading on every indicator is
#'   detected when every indicator also loads on another factor. An explicit
#'   name also supports general factors loading on only some indicators.
#'
#' @details Calculations use fitted loadings, factor
#'   covariances, and the complete fitted residual covariance matrix.
#'   Standardization uses model-implied observed and latent variances.
#'   With `std = "all"`, residual covariances use observed-variable standardization, equivalent to
#'   `cov.std = FALSE` in [lavaan::standardizedSolution()], so tracing products
#'   sum to the model-implied covariance rather than a residual correlation.
#'   Currently supports single-group, single-level CFA models with continuous
#'   indicators and without structural regression or higher-order factors.
#'   Indicators may load on multiple factors; each indicator is drawn once with
#'   an arrow for every fitted loading. Cross-loadings are retained in all
#'   calculations, including the twice-counted cross term in an indicator's
#'   variance when its factors covary. Bifactor layouts do not impose
#'   orthogonality: fitted factor covariances are displayed and used as fitted.
#'   Factor groups sharing indicators are placed together, and coefficient
#'   labels are spaced using their measured size on the active graphics device.
#'   With `std = "lv"`, observed variables and residuals retain their input scales;
#'   only latent factors are rescaled to unit variance. With `std = "estimate"`,
#'   factor scales also remain those of the fit.
#'   In the diagram, double-headed loops beneath observed variables show their
#'   residual variances, not their total variances (these loops appear beside
#'   the indicators in bifactor layouts). Two-headed curves connecting
#'   observed variables show fitted residual covariances. Latent variance loops
#'   appear to the left of the factors. Separate error nodes are omitted;
#'   this changes only the display, not the estimates or calculations.
#'   With `std = "all"`, the requested
#'   covariance between two observed variables is their model-implied
#'   correlation. The tracing key is wrapped and sized to the graphics device.
#' @return Invisibly returns a list containing the plotted coefficients,
#'   factor covariance matrix, residual variances and covariances, joint
#'   model-implied covariance matrix, requested calculations, and individual
#'   tracing contributions.
#'   The `layout`, `nodes`, and `edges` entries describe the drawn arrangement.
#'   Draws the diagram on the active graphics device.
#' @examples
#' \dontrun{
#' #### Trace covariance between two indicators in your fitted model
#' path_diagram(fit, covariance = c("cubes", "word"), tracing = TRUE)
#'
#' #### Specify the same pair in a single comma-separated string
#' path_diagram(fit, covariance = "cubes, word", tracing = TRUE)
#'
#' #### Trace the variance of a latent factor
#' path_diagram(fit, variance = "Visual_LF", tracing = TRUE)
#'
#' #### Choose fully standardized, latent-only standardized, or raw estimates
#' path_diagram(fit, std = "all")
#' path_diagram(fit, std = "lv")
#' path_diagram(fit, std = "estimate")
#'
#' #### Place a fitted general factor below the indicators
#' path_diagram(bifactor_fit, layout = "bifactor", general = "g")
#' }
#' @export

#### Define the function and its settings
path_diagram <- function(
    object,
    digits = 2,
    std = "all",
    variance = NULL,
    covariance = NULL,
    tracing = FALSE,
    layout = c("auto", "tree", "bifactor"),
    general = NULL
) {
  #### Validate the display settings and fitted model
  if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
      !is.finite(digits) || digits < 0 || digits > 15 || digits != floor(digits)) {
    if (is.character(digits) && !is.null(covariance)) {
      stop(paste0("Put both covariance names inside one argument: ",
                  "covariance = c(\"cubes\", \"word\"). ",
                  "An extra unnamed name is treated as `digits`."), call. = FALSE)
    }
    stop("`digits` must be a whole number between 0 and 15.", call. = FALSE)
  }
  if (is.logical(std) && length(std) == 1L && !is.na(std)) {
    std <- if (std) "all" else "estimate"
  }
  if (!is.character(std) || length(std) != 1L || is.na(std) ||
      !std %in% c("all", "lv", "estimate")) {
    stop('`std` must be "all", "lv", or "estimate".', call. = FALSE)
  }
  if (!is.logical(tracing) || length(tracing) != 1L || is.na(tracing)) {
    stop("`tracing` must be TRUE or FALSE.", call. = FALSE)
  }
  layout <- match.arg(layout)
  if (!is.null(general) && (!is.character(general) || length(general) != 1L ||
      is.na(general) || !nzchar(general))) {
    stop("`general` must name one latent factor.", call. = FALSE)
  }
  if (layout == "tree" && !is.null(general)) {
    stop('Use `layout = "auto"` or `"bifactor"` with `general`.', call. = FALSE)
  }

  #### Extract lavaan estimates, including explicitly fitted residual covariances
  is_lavaan <- inherits(object, "lavaan")
  residual_covariance <- NULL
  if (is_lavaan) {
    if (!requireNamespace("lavaan", quietly = TRUE)) {
      stop("Install `lavaan` to plot a lavaan model.", call. = FALSE)
    }
    if (lavaan::lavInspect(object, "ngroups") != 1L) {
      stop("Currently supports single-group lavaan CFA models.", call. = FALSE)
    }
    if (!isTRUE(lavaan::lavInspect(object, "converged"))) {
      stop("The lavaan model must converge before it can be plotted.", call. = FALSE)
    }
    if (length(lavaan::lavNames(object, type = "ov.ord"))) {
      stop("Currently supports continuous lavaan indicators.", call. = FALSE)
    }
    estimates <- lavaan::lavInspect(object, "est")
    parameters <- lavaan::parameterTable(object)
    if (is.null(estimates$lambda) || is.null(estimates$psi) ||
        is.null(estimates$theta) || any(parameters$op == "~") ||
        (!is.null(estimates$beta) && any(abs(estimates$beta) > 1e-10))) {
      stop(paste0("Currently supports single-level CFA models without ",
                  "structural regression or higher-order factors."), call. = FALSE)
    }
    residual_covariance <- as.matrix(estimates$theta)
    model_covariance <- estimates$lambda %*% estimates$psi %*% t(estimates$lambda) +
      residual_covariance
    fit <- list(Pattern = estimates$lambda, loadings = estimates$lambda,
                Phi = estimates$psi, r = model_covariance)
  } else {
    #### Retain the original input extraction for existing fitted objects
    if (!inherits(object, "CFA") || !is.list(object$cfa) ||
        is.null(object$cfa$Pattern) || is.null(object$cfa$r)) {
      stop("`object` must be returned by lavaan::cfa().", call. = FALSE)
    }
    fit <- object$cfa
  }

  #### Extract the original fitted pattern coefficients and input matrix
  loading <- as.matrix(fit$Pattern)

  #### Preserve loading signs in single-factor structure loadings
  if (ncol(loading) == 1L && !is.null(fit$loadings)) {
    loading <- as.matrix(fit$loadings)
  }
  input <- as.matrix(fit$r)
  p <- nrow(loading)
  k <- ncol(loading)
  observed <- rownames(loading)
  latent <- colnames(loading)

  #### Recover single-factor names retained on the structure loading matrix
  if (is.null(observed)) observed <- rownames(fit$loadings)
  if (is.null(latent)) latent <- colnames(fit$loadings)
  dimnames(loading) <- list(observed, latent)
  if (!is.numeric(loading) || !length(loading) || any(!is.finite(loading)) ||
      !is.numeric(input) || !identical(dim(input), c(p, p)) ||
      any(!is.finite(input))) {
    stop("The model must contain finite numeric pattern and input matrices.", call. = FALSE)
  }
  if (is.null(observed) || is.null(latent) ||
      anyNA(c(observed, latent)) || any(!nzchar(c(observed, latent))) ||
      anyDuplicated(c(observed, latent))) {
    stop("Observed and latent variables must have distinct, nonempty names.", call. = FALSE)
  }
  if (!is.null(general) && !general %in% latent) {
    stop("`general` must name one latent factor. Available factors: ",
         paste(latent, collapse = ", "), call. = FALSE)
  }

  #### Preserve model topology, including freely estimated loadings equal to zero
  loading_paths <- abs(loading) > 1e-10
  if (is_lavaan) {
    loading_paths[,] <- FALSE
    rows <- which(parameters$op == "=~" &
                    (parameters$free > 0 | abs(parameters$est) > 1e-10))
    for (row in rows) {
      loading_paths[parameters$rhs[row], parameters$lhs[row]] <- TRUE
    }
  }
  if (is.null(general) && layout != "tree" && k > 1L) {
    candidates <- which(colSums(loading_paths) == p)
    candidates <- candidates[vapply(candidates, function(i)
      all(rowSums(loading_paths[, -i, drop = FALSE]) > 0), logical(1))]
    if (length(candidates) == 1L) general <- latent[candidates]
  }
  if (layout == "bifactor" && is.null(general)) {
    stop(paste0("Could not identify a unique general factor. ",
                "Set `general` to its latent variable name."), call. = FALSE)
  }
  if (!is.null(general) && (k < 2L || !any(loading_paths[, general]))) {
    stop("A bifactor layout requires a general factor with indicators and at least one other factor.",
         call. = FALSE)
  }
  layout <- if (is.null(general)) "tree" else "bifactor"
  if (!is.null(rownames(input)) && !is.null(colnames(input))) {
    if (!all(observed %in% rownames(input)) || !all(observed %in% colnames(input))) {
      stop("The input matrix names do not match the observed variables.", call. = FALSE)
    }
    input <- input[observed, observed, drop = FALSE]
  }

  #### Extract factor covariances, including the single-factor case
  phi <- fit$Phi
  if (is.null(phi)) {
    phi <- diag(k)
  }
  phi <- as.matrix(phi)
  if (!is.numeric(phi) || !identical(dim(phi), c(k, k)) ||
      any(!is.finite(phi))) {
    stop("The model contains an invalid factor covariance matrix.", call. = FALSE)
  }
  if (!is.null(rownames(phi)) && !is.null(colnames(phi))) {
    if (!all(latent %in% rownames(phi)) || !all(latent %in% colnames(phi))) {
      stop("The factor covariance names do not match the factors.", call. = FALSE)
    }
    phi <- phi[latent, latent, drop = FALSE]
  }
  dimnames(phi) <- list(latent, latent)
  if (!isTRUE(all.equal(input, t(input), check.attributes = FALSE)) ||
      !isTRUE(all.equal(phi, t(phi), check.attributes = FALSE)) ||
      any(diag(input) <= 0) || any(diag(phi) <= 0) ||
      min(eigen(phi, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
    stop("Input and factor matrices must be symmetric with valid variances.", call. = FALSE)
  }

  #### Require retained observed-variable scales for unstandardized estimates
  if (!is_lavaan && std == "estimate" && all(abs(diag(input) - 1) < 1e-8)) {
    stop(paste0(
      '`std = "estimate"` requires retained observed-variable scales. ',
      "Fit using a covariance matrix or retain the original observed-variable scales. ",
      "A unit-diagonal matrix is treated as correlation-only."
    ), call. = FALSE)
  }

  #### Retain fitted residual covariances or derive unique residual variances
  common <- loading %*% phi %*% t(loading)
  if (is.null(residual_covariance)) {
    residual_covariance <- diag(diag(input) - diag(common), nrow = p)
  }
  dimnames(residual_covariance) <- list(observed, observed)
  residual <- diag(residual_covariance)
  tolerance <- 1e-8 * max(diag(input))
  if (any(residual < -tolerance) ||
      min(eigen(residual_covariance, symmetric = TRUE, only.values = TRUE)$values) < -tolerance) {
    stop("The fitted model contains invalid residual variances or covariances.", call. = FALSE)
  }
  diag(residual_covariance) <- pmax(residual, 0)

  #### Standardize latent factors, and optionally the observed variables
  if (std != "estimate") {
    observed_sd <- if (std == "all") sqrt(diag(input)) else rep(1, p)
    latent_sd <- sqrt(diag(phi))
    loading <- sweep(sweep(loading, 1, observed_sd, "/"), 2, latent_sd, "*")
    phi <- phi / outer(latent_sd, latent_sd)
    residual_covariance <- residual_covariance / outer(observed_sd, observed_sd)
  }
  residual <- diag(residual_covariance)
  names(residual) <- observed

  #### Remove numerical noise so zero paths do not appear in the tracing key
  loading[abs(loading) <= 1e-10] <- 0
  phi[abs(phi) <= 1e-10] <- 0
  residual_covariance[abs(residual_covariance) <= 1e-10] <- 0

  #### Construct the joint model-implied covariance matrix
  variable_names <- c(latent, observed)
  weights <- rbind(diag(k), loading)
  dimnames(weights) <- list(variable_names, latent)
  implied <- weights %*% phi %*% t(weights)
  observed_indices <- k + seq_len(p)
  implied[observed_indices, observed_indices] <-
    implied[observed_indices, observed_indices, drop = FALSE] + residual_covariance
  dimnames(implied) <- list(variable_names, variable_names)

  #### Validate named requests and allow a comma-separated covariance pair
  if (is.character(covariance) && length(covariance) == 1L &&
      !is.na(covariance) && grepl(",", covariance, fixed = TRUE)) {
    covariance <- trimws(strsplit(covariance, ",", fixed = TRUE)[[1]])
  }
  for (argument in c("variance", "covariance")) {
    request <- get(argument)
    required_length <- if (argument == "variance") 1L else 2L
    if (!is.null(request) && (!is.character(request) ||
        length(request) != required_length || anyNA(request) ||
        !all(request %in% variable_names))) {
      example <- if (argument == "covariance")
        ' Use covariance = c("cubes", "word") or covariance = "cubes, word".' else
        ' Use variance = "Visual_LF".'
      stop("`", argument, "` must name ", required_length,
           " variable(s) in the model.", example, " Available names: ",
           paste(variable_names, collapse = ", "), call. = FALSE)
    }
  }
  if (tracing && is.null(variance) && is.null(covariance)) {
    stop("Choose `variance` or `covariance` when `tracing = TRUE`.", call. = FALSE)
  }

  #### Create requests for model-implied variance and covariance calculations
  requests <- list()
  if (!is.null(variance)) requests$variance <- c(variance, variance)
  if (!is.null(covariance)) requests$covariance <- covariance
  calculations <- data.frame(
    calculation = character(), from = character(), to = character(),
    value = numeric(), stringsAsFactors = FALSE
  )
  traces <- data.frame(
    calculation = character(), path = character(), expression = character(),
    value = numeric(), stringsAsFactors = FALSE
  )
  trace_edges <- list()

  #### Identify the diagram edge for each loading or factor relationship
  loading_edge <- function(variable, factor) {
    if (variable %in% latent) return(character())
    paste0("L", match(variable, observed), ":", factor)
  }
  factor_edge <- function(i, j) {
    if (i == j) return(paste0("V", i))
    paste0("C", min(i, j), ":", max(i, j))
  }
  number <- function(x) formatC(x, format = "f", digits = digits)

  #### Enumerate valid paths and their products without rounding calculations
  for (calculation in names(requests)) {
    pair <- requests[[calculation]]
    a <- pair[1]
    b <- pair[2]
    calculations[nrow(calculations) + 1L, ] <- list(
      calculation, a, b, unname(implied[a, b])
    )
    for (i in seq_len(k)) {
      for (j in seq_len(k)) {
        if (weights[a, i] == 0 || weights[b, j] == 0 || phi[i, j] == 0 ||
            (a == b && j < i)) next
        multiplier <- if (a == b && i != j) 2 else 1
        factors <- c(weights[a, i], phi[i, j], weights[b, j])
        path <- c(a, latent[i], latent[j], b)
        path <- path[c(TRUE, path[-1] != path[-length(path)])]
        traces[nrow(traces) + 1L, ] <- list(
          calculation, paste(path, collapse = " -> "),
          paste0(if (multiplier == 2) "2 * " else "",
                 paste(number(factors), collapse = " * ")),
          multiplier * prod(factors)
        )
        trace_edges[[nrow(traces)]] <- unique(c(
          loading_edge(a, i), factor_edge(i, j), loading_edge(b, j)
        ))
      }
    }

    #### Include the unique residual contribution to an observed variance
    if (a == b && a %in% observed && residual[a] != 0) {
      index <- match(a, observed)
      traces[nrow(traces) + 1L, ] <- list(
        calculation, paste0("Residual variance loop at ", a),
        number(residual[a]), unname(residual[a])
      )
      trace_edges[[nrow(traces)]] <- paste0("U", index)
    }

    #### Include the fitted residual covariance between two observed variables
    if (a != b && all(c(a, b) %in% observed) && residual_covariance[a, b] != 0) {
      i <- match(a, observed)
      j <- match(b, observed)
      traces[nrow(traces) + 1L, ] <- list(
        calculation, paste0(a, " <-> ", b, " (residual covariance)"),
        number(residual_covariance[a, b]),
        unname(residual_covariance[a, b])
      )
      trace_edges[[nrow(traces)]] <- paste0("R", min(i, j), ":", max(i, j))
    }
  }

  #### Keep connected factor groups together before placing shared indicators
  general_index <- match(general, latent)
  upper <- setdiff(seq_len(k), general_index)
  shared <- crossprod(loading_paths[, upper, drop = FALSE])
  diag(shared) <- 0
  remaining <- seq_along(upper)
  factor_order <- integer()
  while (length(remaining)) {
    queue <- remaining[1L]
    while (length(queue)) {
      next_factor <- queue[1L]
      queue <- queue[-1L]
      if (!next_factor %in% remaining) next
      factor_order <- c(factor_order, next_factor)
      remaining <- setdiff(remaining, next_factor)
      neighbors <- remaining[shared[next_factor, remaining] > 0]
      neighbors <- neighbors[order(-shared[next_factor, neighbors], neighbors)]
      queue <- unique(c(queue, neighbors))
    }
  }
  upper <- upper[factor_order]
  #### Membership determines positions, so changing estimate scales cannot move nodes
  membership <- loading_paths[, upper, drop = FALSE]
  group <- as.vector(membership %*% seq_along(upper)) / pmax(rowSums(membership), 1)
  group[rowSums(membership) == 0] <- (length(upper) + 1) / 2
  order_observed <- order(group, seq_len(p))
  observed_x <- numeric(p)
  observed_x[order_observed] <- (seq_len(p) - 0.5) / p
  factor_x <- factor_y <- numeric(k)
  factor_x[upper] <- (seq_along(upper) - 0.5) / length(upper)
  factor_y[upper] <- 0.95
  if (length(general_index)) {
    factor_x[general_index] <- mean(observed_x[loading_paths[, general_index]])
    factor_y[general_index] <- 0.07
  }
  observed_y <- if (layout == "bifactor") 0.52 else 0.30
  box_bottom <- observed_y - 0.045
  box_top <- observed_y + 0.045
  observed_width <- min(0.12, 0.8 / p)
  factor_width <- min(0.18, 0.8 / length(upper))
  nodes <- data.frame(name = variable_names, type = c(rep("latent", k), rep("observed", p)),
                      x = c(factor_x, observed_x), y = c(factor_y, rep(observed_y, p)),
                      width = c(rep(factor_width, k), rep(observed_width, p)),
                      height = c(rep(0.08, k), rep(0.09, p)), stringsAsFactors = FALSE)
  edges <- list()

  #### Store curved paths with arrow directions and coefficient label positions
  add_edge <- function(id, x, y, value, both = FALSE, type = "path",
                       from = NULL, to = NULL, label_position = 0.5) {
    edges[[id]] <<- list(x = x, y = y, value = value, both = both, type = type,
                        from = from, to = to, label_position = label_position)
  }
  t <- seq(0, 1, length.out = 80)
  for (i in seq_len(k)) {
    for (j in seq_len(p)) {
      if (!loading_paths[j, i]) next
      below <- i %in% general_index
      incoming <- upper[loading_paths[j, upper]]
      attachment <- if (!below && length(incoming) > 1L)
        (match(i, incoming) - (length(incoming) + 1) / 2) *
          observed_width * 0.65 / (length(incoming) - 1) else 0
      end_x <- observed_x[j] + attachment
      end_y <- if (below) box_bottom else box_top
      dx <- end_x - factor_x[i]
      dy <- end_y - factor_y[i]
      boundary <- 1 / sqrt((dx / (factor_width / 2))^2 + (dy / 0.04)^2)
      add_edge(paste0("L", j, ":", i),
               c(factor_x[i] + boundary * dx, end_x),
               c(factor_y[i] + boundary * dy, end_y), loading[j, i],
               type = "loading", from = latent[i], to = observed[j],
               label_position = if (below) 0.66 else if (length(incoming) > 1L)
                 0.30 + 0.18 * match(i, upper) / length(upper) else 0.5)
    }
    #### Attach a double-headed variance loop to the left side of each latent factor
    theta <- seq(-0.65 * pi, 0.65 * pi, length.out = 80)
    loop_width <- min(0.03, 0.15 / k)
    loop_height <- 0.03
    attachment_y <- loop_height * sin(0.65 * pi)
    attachment_x <- factor_width / 2 * sqrt(1 - (attachment_y / 0.04)^2)
    loop_center <- factor_x[i] - attachment_x + loop_width * cos(0.65 * pi)
    add_edge(paste0("V", i),
             loop_center - loop_width * cos(theta),
             factor_y[i] + loop_height * sin(theta), phi[i, i], TRUE,
             type = "latent_variance", from = latent[i], to = latent[i])
    if (i < k) {
      for (j in seq.int(i + 1L, k)) {
        if (abs(phi[i, j]) <= 1e-10) next
        if (any(c(i, j) %in% general_index)) {
          #### Route general/specific covariances outside the indicator row
          specific <- setdiff(c(i, j), general_index)
          side <- if (factor_x[specific] <= 0.5) -1 else 1
          start_x <- factor_x[specific] + side * factor_width / 2
          end_x <- factor_x[general_index] + side * factor_width / 2
          lane <- if (side < 0) -0.16 - 0.12 * match(specific, upper) else
            1.16 + 0.12 * match(specific, upper)
          outside <- (lane - 0.125 * (start_x + end_x)) / 0.75
          x <- (1 - t)^3 * start_x + 3 * (1 - t)^2 * t * outside +
            3 * (1 - t) * t^2 * outside + t^3 * end_x
          y <- (1 - t)^3 * 0.95 + 3 * (1 - t)^2 * t * 0.65 +
            3 * (1 - t) * t^2 * 0.07 + t^3 * 0.07
        } else {
          x <- (1 - t) * factor_x[i] + t * factor_x[j]
          y <- 0.99 + (0.13 + 0.12 * abs(match(j, upper) - match(i, upper)) /
                        length(upper)) * 4 * t * (1 - t)
        }
        add_edge(paste0("C", i, ":", j), x, y, phi[i, j], TRUE,
                 type = "factor_covariance", from = latent[i], to = latent[j])
      }
    }
  }

  #### Keep residual loops clear of general-factor arrows from below
  for (j in seq_len(p)) {
    theta <- seq(-0.65 * pi, 0.65 * pi, length.out = 80)
    loop_width <- min(0.025, observed_width * 0.3)
    if (layout == "bifactor") {
      loop_width <- min(0.016, (1 / p - observed_width) / 3)
      loop_center <- observed_x[j] + observed_width / 2 - loop_width * cos(0.65 * pi)
      x <- loop_center + loop_width * cos(theta)
      y <- observed_y + 0.025 * sin(theta)
    } else {
      loop_center <- box_bottom + 0.055 * cos(0.65 * pi)
      x <- observed_x[j] + loop_width * sin(theta)
      y <- loop_center - 0.055 * cos(theta)
    }
    add_edge(paste0("U", j), x, y, residual[j], TRUE,
             type = "residual_variance", from = observed[j], to = observed[j])
  }

  #### Draw residual covariance curves below the boxes and their variance loops
  if (p > 1L) {
    for (i in seq_len(p - 1L)) {
      for (j in seq.int(i + 1L, p)) {
        if (abs(residual_covariance[i, j]) <= 1e-10) next
        add_edge(paste0("R", i, ":", j),
                 (1 - t) * observed_x[i] + t * observed_x[j],
                 box_bottom - (if (layout == "bifactor") 0.10 else 0.18 +
                                 0.04 * (j - i) / p) * 4 * t * (1 - t),
                 residual_covariance[i, j], TRUE, type = "residual_covariance",
                 from = observed[i], to = observed[j])
      }
    }
  }

  #### Initialize the graphics device before measuring text in physical units
  previous <- graphics::par(mar = c(0.25, 0.3, 1.8, 0.3))
  on.exit(graphics::par(previous), add = TRUE)
  graphics::plot.new()
  plot_inches <- graphics::par("pin")
  #### Set number size, word size, font weight, and text shade
  coefficient_cex <- min(0.85, max(0.65, plot_inches[1] / 9))
  node_cex <- min(0.9, max(0.65, plot_inches[1] / 8))
  number_font <- 2
  word_font <- 1
  text_color <- "black"

  #### Assemble calculation totals, tracing explanations, and the residual guide
  footer <- list()
  add_footer <- function(text, color = text_color, font = word_font, swatch = FALSE, scale = 1) {
    footer[[length(footer) + 1L]] <<- list(
      text = text, color = color, font = font, swatch = swatch, scale = scale
    )
  }
  for (i in seq_len(nrow(calculations))) {
    item <- calculations[i, ]
    label <- if (item$calculation == "variance") paste0("Var(", item$from, ")") else
      paste0("Cov(", item$from, ", ", item$to, ")")
    add_footer(paste0(label, " = ", number(item$value)), font = number_font)
  }
  colors <- if (nrow(traces)) grDevices::hcl.colors(nrow(traces), "Dark 3") else character()
  if (tracing && nrow(traces)) {
    add_footer("Tracing key: multiply along each path, then add the contributions.", font = word_font)
    for (i in seq_len(nrow(traces))) {
      path <- traces$path[i]
      if (path %in% latent) path <- paste0("Variance loop at ", path)
      prefix <- if (nrow(calculations) > 1L) paste0(traces$calculation[i], ": ") else ""
      add_footer(paste0(prefix, path, " | ", traces$expression[i], " = ", number(traces$value[i])),
                 color = colors[i], font = word_font, swatch = TRUE)
    }
  } else if (tracing) {
    add_footer("No contributing paths: the requested covariance is zero.")
  }
  add_footer(paste0("Double-headed loops ",
                    if (layout == "bifactor") "beside" else "beneath",
                    " observed variables show residual variances."),
             color = "#555555", scale = 0.9)
  if (any(residual_covariance[row(residual_covariance) != col(residual_covariance)] != 0)) {
    add_footer("Two-headed curves between observed variables show residual covariances.",
               color = "#555555", scale = 0.9)
  }

  #### Wrap text to measured device width, including unusually long variable names
  wrap_text <- function(text, width, cex, font) {
    lines <- character()
    remaining <- text
    while (graphics::strwidth(remaining, units = "inches", cex = cex, font = font) > width) {
      indices <- seq_len(nchar(remaining))
      fits <- indices[graphics::strwidth(substring(remaining, 1, indices), units = "inches",
                                        cex = cex, font = font) <= width]
      cut <- if (length(fits)) max(fits) else 1L
      spaces <- gregexpr(" ", substr(remaining, 1, cut), fixed = TRUE)[[1]]
      spaces <- spaces[spaces > 0]
      if (length(spaces) && max(spaces) > cut / 2) cut <- max(spaces)
      lines <- c(lines, trimws(substr(remaining, 1, cut)))
      remaining <- trimws(substring(remaining, cut + 1L))
    }
    c(lines, remaining)
  }

  #### Fit every key line inside a footer whose height is measured in inches
  footer_cex <- 0.8
  for (attempt in seq_len(40)) {
    footer_lines <- list()
    for (entry in footer) {
      available <- max(0.1, plot_inches[1] - 0.08 - if (entry$swatch) 0.26 else 0)
      lines <- wrap_text(entry$text, available, footer_cex * entry$scale, entry$font)
      for (j in seq_along(lines)) {
        item <- entry
        item$text <- lines[j]
        item$first <- j == 1L
        footer_lines[[length(footer_lines) + 1L]] <- item
      }
    }
    line_inches <- graphics::strheight("Mg", units = "inches", cex = footer_cex) * 1.65
    footer_inches <- (length(footer_lines) + 0.9) * line_inches + 0.06
    if (footer_inches <= plot_inches[2] * 0.45) break
    footer_cex <- footer_cex * 0.92
  }
  footer_fraction <- footer_inches / plot_inches[2]
  lower <- -1.3 * footer_fraction / (1 - footer_fraction)
  #### Reserve horizontal space for side variance loops and their bold coefficients
  loop_left <- min(vapply(edges[paste0("V", seq_len(k))], function(edge) min(edge$x), numeric(1)))
  edge_x <- unlist(lapply(edges, function(edge) range(edge$x)), use.names = FALSE)
  graphics::plot.window(xlim = c(min(-0.03, loop_left - 0.085, edge_x - 0.03),
                                max(1.03, edge_x + 0.03)),
                         ylim = c(lower, 1.3), xaxs = "i", yaxs = "i")
  scale_label <- switch(std, all = "standardized: Std.all",
                        lv = "latent standardized: Std.lv", estimate = "unstandardized: Estimate")
  graphics::title(main = paste0("CFA path diagram (", scale_label, ")"),
                   cex.main = min(1.2, max(0.8, plot_inches[1] / 5.5)),
                   font.main = word_font, col.main = text_color)

  #### Fit coefficient labels between nodes and previously placed labels
  bounds <- function(x, y, width, height) c(x - width / 2, x + width / 2,
                                           y - height / 2, y + height / 2)
  node_bounds <- t(vapply(seq_len(nrow(nodes)), function(i)
    bounds(nodes$x[i], nodes$y[i], nodes$width[i], nodes$height[i]), numeric(4)))
  label_bounds <- matrix(numeric(), ncol = 4)
  overlap_count <- function(box, others) {
    sum(box[1] < others[, 2] & box[2] > others[, 1] &
          box[3] < others[, 4] & box[4] > others[, 3])
  }
  edge_types <- vapply(edges, `[[`, "", "type")
  loops <- edge_types %in% c("latent_variance", "residual_variance")
  for (id in names(edges)[c(which(loops), which(!loops))]) {
    edge <- edges[[id]]
    text <- number(edge$value)
    fractions <- if (edge$type %in% c("loading", "factor_covariance", "residual_covariance"))
      unique(pmax(0.2, pmin(0.8, edge$label_position +
                              c(0, .1, -.1, .2, -.2, .3, -.3)))) else 0.5
    best <- NULL
    best_penalty <- Inf
    for (label_cex in unique(pmax(0.5, coefficient_cex * c(1, .92, .84, .76)))) {
      width <- graphics::strwidth(text, cex = label_cex, font = number_font) * 1.15
      height <- graphics::strheight(text, cex = label_cex, font = number_font) * 1.25
      for (fraction in fractions) {
        if (edge$type == "latent_variance") {
          x <- min(edge$x) - graphics::xinch(0.025) - width / 2
          y <- mean(range(edge$y))
        } else if (edge$type == "residual_variance") {
          x <- mean(range(edge$x))
          y <- if (layout == "bifactor") max(box_top, edge$y) + height else
            min(edge$y) - height
        } else {
          positions <- seq(0, 1, length.out = length(edge$x))
          x <- stats::approx(positions, edge$x, xout = fraction)$y
          y <- stats::approx(positions, edge$y, xout = fraction)$y + height * 0.8
        }
        box <- bounds(x, y, width, height)
        padded <- box + c(-1, 1, -1, 1) *
          c(rep(graphics::xinch(0.02), 2), rep(graphics::yinch(0.015), 2))
        usr <- graphics::par("usr")
        penalty <- overlap_count(padded, rbind(node_bounds, label_bounds)) +
          sum(box[c(1, 3)] < usr[c(1, 3)] | box[c(2, 4)] > usr[c(2, 4)])
        if (penalty < best_penalty) {
          best <- list(x = x, y = y, width = width, height = height,
                       cex = label_cex, text = text)
          best_penalty <- penalty
        }
        if (penalty == 0) break
      }
      if (best_penalty == 0) break
    }
    edges[[id]]$label <- best
    label_bounds <- rbind(label_bounds, bounds(best$x, best$y, best$width, best$height))
  }

  #### Draw paths with arrows and their coefficient labels
  draw_edge <- function(edge, color, width = 1, offset = 0, label = FALSE) {
    x <- edge$x + offset
    y <- edge$y
    n <- length(x)

    #### Label residual variances below their loops and factor variances to the left
    if (label) {
      position <- edge$label
      box <- bounds(position$x, position$y, position$width, position$height)
      graphics::rect(box[1], box[3], box[2], box[4], col = "white", border = NA)
      graphics::text(position$x, position$y, position$text, cex = position$cex,
                     font = number_font, col = text_color)
      return(invisible(NULL))
    }
    graphics::lines(x, y, col = color, lwd = width)
    arrow_inches <- min(0.06, plot_inches[1] / 120)
    graphics::arrows(x[n - 1], y[n - 1], x[n], y[n], length = arrow_inches, col = color, lwd = width)
    if (edge$both) {
      graphics::arrows(x[2], y[2], x[1], y[1], length = arrow_inches, col = color, lwd = width)
    }
  }
  for (edge in edges) draw_edge(edge, if (tracing) "#A0A5AD" else "#68717D", width = 1.1)

  #### Overlay a distinct color for each contribution and retain shared paths
  if (tracing) {
    for (i in seq_len(nrow(traces))) {
      for (id in trace_edges[[i]]) {
        if (is.null(edges[[id]])) next
        draw_edge(edges[[id]], colors[i], 2.5,
                  offset = (i - (nrow(traces) + 1) / 2) * min(0.002, 0.018 / nrow(traces)))
      }
    }
  }
  for (edge in edges) draw_edge(edge, "#333333", label = TRUE)

  #### Draw latent ellipses and observed boxes over the paths
  ellipse <- function(x, y, width, height, label, fill) {
    angle <- seq(0, 2 * pi, length.out = 80)
    graphics::polygon(x + width / 2 * cos(angle), y + height / 2 * sin(angle),
                      col = fill, border = "#333333")
    label_cex <- min(node_cex, node_cex * width * 0.85 /
                      max(graphics::strwidth(label, cex = node_cex, font = word_font), 0.001))
    graphics::text(x, y, label, cex = label_cex, font = word_font, col = text_color)
  }
  for (i in seq_len(k)) ellipse(factor_x[i], factor_y[i], factor_width, 0.08, latent[i], "#DDEAF5")
  for (j in seq_len(p)) {
    graphics::rect(observed_x[j] - observed_width / 2, box_bottom,
                   observed_x[j] + observed_width / 2, box_top,
                   col = "#FFFFFF", border = "#333333")
    label_cex <- min(node_cex, node_cex * observed_width * 0.9 /
                      max(graphics::strwidth(observed[j], cex = node_cex, font = word_font), 0.001))
    graphics::text(observed_x[j], observed_y, observed[j], cex = label_cex, font = word_font, col = text_color)
  }

  #### Place each wrapped key line inside the measured footer without clipping
  footer_x <- graphics::par("usr")[1] + graphics::xinch(0.04)
  footer_y <- -graphics::yinch(0.04)
  for (i in seq_along(footer_lines)) {
    item <- footer_lines[[i]]
    y <- footer_y - graphics::yinch((i - 0.5) * line_inches)
    indent <- if (item$swatch) graphics::xinch(0.26) else 0
    if (item$swatch && item$first) {
      graphics::segments(footer_x, y, footer_x + graphics::xinch(0.19), y,
                         col = item$color, lwd = 2.5)
    }
    graphics::text(footer_x + indent, y, item$text, adj = c(0, 0.5),
                   cex = footer_cex * item$scale, font = item$font,
                   col = text_color)
  }

  #### Return full-precision estimates and tracing contributions invisibly
  traces$color <- colors
  traces$edges <- I(trace_edges)
  invisible(list(
    std = std, loadings = loading, factor_covariance = phi,
    residual_variance = residual, residual_covariance = residual_covariance,
    implied_covariance = implied,
    calculations = calculations, traces = traces,
    layout = layout, nodes = nodes, edges = edges
  ))
}
