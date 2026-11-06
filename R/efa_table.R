#' Create a shaded EFA loading table
#'
#' @param object A fitted model containing a numeric loading matrix.
#' @param shadestyle Either "posneg" for signed shading (meaning two
#'   unique colors for positive and negative) or "abs"
#'   for shading by absolute magnitude.
#' @param digits Number of decimal places to display.
#' @param negcol Color for negative loadings.
#' @param zerocol Color for zero loadings.
#' @param poscol Color for positive loadings.
#'
#' @return Outputs a gt table.
#' @export

#### Define the function and its settings
efa_table <- function(
    object,
    shadestyle = "posneg",
    digits = 2,
    negcol = "#6BAED6",
    zerocol = "#FFFFFF",
    poscol = "#F28E2B"
) {

  #### Validate the requested shading style
  shadestyle <- match.arg(
    shadestyle,
    choices = c("posneg", "abs")
  )

  #### Extract loadings from the supplied model
  Load <- unclass(object$loadings)

  #### Calculate the shading limit
  limit <- max(abs(Load), na.rm = TRUE)

  if (!is.finite(limit) || limit == 0) {
    limit <- 1
  }

  #### Define each shading style using the supplied colors
  if (shadestyle == "abs") {
    shade_palette <- c(zerocol, poscol)
    color_range <- c(0, limit)

  } else if (shadestyle == "posneg") {
    shade_palette <- c(negcol, zerocol, poscol)
    color_range <- c(-limit, limit)

  } else {
    stop("No palette has been defined for this shadestyle.")
  }

  #### Translate loading values into colors
  loading_colors <- function(x) {
    if (shadestyle == "abs") {
      x <- abs(x)
    }

    scales::col_numeric(
      palette = shade_palette,
      domain = color_range,
      na.color = "#EEEEEE"
    )(x)
  }

  #### Build and return the table
  gt::gt(
    as.data.frame(Load),
    rownames_to_stub = TRUE
  ) |>
    gt::fmt_number(decimals = digits) |>
    gt::data_color(
      fn = loading_colors,
      autocolor_text = FALSE
    )
}
