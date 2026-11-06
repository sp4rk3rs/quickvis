#' Compare EFA loading tables side by side
#'
#' @param ... Fitted models, optionally named to supply model headings.
#' @param model_gap Space between models, in pixels.
#' @inheritParams efa_table
#'
#' @return A gt table with model headings and a shared color scale.
#' @export

#### Define a function for comparing multiple EFA models
compare_efa <- function(
    ...,
    shadestyle = "posneg",
    digits = 2,
    model_gap = 24,
    negcol = "#6BAED6",
    zerocol = "#FFFFFF",
    poscol = "#F28E2B"
) {

  #### Collect the supplied models
  models <- list(...)

  if (length(models) == 0) {
    stop("Supply at least one factor analysis result.")
  }

  #### Create titles for unnamed models
  titles <- names(models)

  if (is.null(titles)) {
    titles <- rep("", length(models))
  }

  unnamed <- is.na(titles) | !nzchar(titles)
  titles[unnamed] <- paste("Model", which(unnamed))

  #### Extract each model's loadings
  matrices <- lapply(models, function(model) {
    Load <- unclass(model$loadings)

    if (!is.matrix(Load) || !is.numeric(Load)) {
      stop("Each model must contain a numeric loading matrix.")
    }

    Load
  })

  #### Use the first model's item order
  items <- rownames(matrices[[1]])

  if (
    is.null(items) ||
    anyNA(items) ||
    any(!nzchar(items)) ||
    anyDuplicated(items) > 0
  ) {
    stop("The first model must have unique, nonempty item names.")
  }

  #### Match every model's rows to those item names
  for (i in seq_along(matrices)) {
    model_items <- rownames(matrices[[i]])

    if (
      is.null(model_items) ||
      anyDuplicated(model_items) > 0 ||
      !setequal(items, model_items)
    ) {
      stop(
        paste0(
          "'", titles[i],
          "' must contain the same unique item names as the first model."
        )
      )
    }

    matrices[[i]] <- matrices[[i]][items, , drop = FALSE]
  }

  #### Combine the loading matrices side by side
  combined <- do.call(cbind, matrices)

  #### Preserve factor labels while making internal names unique
  factor_labels <- colnames(combined)
  column_ids <- paste0("loading_", seq_len(ncol(combined)))
  colnames(combined) <- column_ids

  #### Reuse efa_table with one shared shading scale
  result <- efa_table(
    object = list(loadings = combined),
    shadestyle = shadestyle,
    digits = digits,
    negcol = negcol,
    zerocol = zerocol,
    poscol = poscol
  )

  #### Restore the displayed factor labels
  result <- gt::cols_label(
    result,
    .list = stats::setNames(
      as.list(factor_labels),
      column_ids
    )
  )

  #### Add a heading above each model's columns
  position <- 0

  for (i in seq_along(matrices)) {
    positions <- position + seq_len(ncol(matrices[[i]]))

    result <- gt::tab_spanner(
      result,
      label = titles[i],
      columns = tidyselect::all_of(column_ids[positions]),
      id = paste0("model_", i)
    )

    position <- position + ncol(matrices[[i]])
  }

  #### Add white space between models
  if (length(matrices) > 1 && model_gap > 0) {

    #### Find the first factor column of each subsequent model
    model_ends <- cumsum(
      vapply(matrices, ncol, integer(1))
    )

    gap_columns <- column_ids[
      model_ends[-length(model_ends)] + 1
    ]

    gap_spanners <- paste0(
      "model_",
      seq.int(2, length(matrices))
    )

    #### Apply white separators to the body and headings
    result <- gt::tab_style(
      result,
      style = gt::cell_borders(
        sides = "left",
        color = "#FFFFFF",
        weight = gt::px(model_gap)
      ),
      locations = list(
        gt::cells_body(
          columns = tidyselect::all_of(gap_columns)
        ),
        gt::cells_column_labels(
          columns = tidyselect::all_of(gap_columns)
        ),
        gt::cells_column_spanners(
          spanners = tidyselect::all_of(gap_spanners)
        )
      )
    )
  }

  #### Return the comparison table
  result
}
