#### Reproduce the shared Visual/Speed indicator shown in the RStudio plot pane
library(quickvis)
if (requireNamespace("lavaan", quietly = TRUE)) {
  data <- lavaan::HolzingerSwineford1939
  observed <- c("visual", "cubes", "lozenge", "parag", "sentence", "word",
                "spd_add", "spd_cnt", "spd_disc")
  names(data)[match(paste0("x", 1:9), names(data))] <- observed
  model <- 'Visual_LF =~ NA*visual + cubes + lozenge + spd_disc
            Verbal =~ NA*parag + sentence + word
            Speed =~ NA*spd_add + spd_cnt + spd_disc
            Visual_LF ~~ 1*Visual_LF
            Verbal ~~ 1*Verbal
            Speed ~~ 1*Speed'
  fit <- lavaan::cfa(model, data = data, missing = "ML", estimator = "MLR")
  intersects <- function(a, b) {
    a[1] < b[2] && a[2] > b[1] && a[3] < b[4] && a[4] > b[3]
  }
  label_box <- function(label) {
    c(label$x - label$width / 2, label$x + label$width / 2,
      label$y - label$height / 2, label$y + label$height / 2)
  }
  for (size in list(c(600, 380), c(1200, 800))) {
    grDevices::png(tempfile(fileext = ".png"), width = size[1], height = size[2], res = 96)
    for (scale in c("lv", "all")) {
      output <- path_diagram(fit, std = scale, digits = 3,
                             covariance = c("word", "cubes"), tracing = TRUE)
      latent_nodes <- output$nodes[output$nodes$type == "latent", ]
      item_nodes <- output$nodes[output$nodes$type == "observed", ]
      stopifnot(identical(latent_nodes$name[order(latent_nodes$x)],
                          c("Visual_LF", "Speed", "Verbal")),
                item_nodes$x[item_nodes$name == "spd_disc"] >
                  max(item_nodes$x[item_nodes$name %in% c("visual", "cubes", "lozenge")]),
                item_nodes$x[item_nodes$name == "spd_disc"] <
                  min(item_nodes$x[item_nodes$name %in% c("spd_add", "spd_cnt")]),
                abs(sum(output$traces$value) - output$calculations$value) < 1e-8)
      boxes <- lapply(output$edges, function(edge) label_box(edge$label))
      for (i in seq_along(boxes)) {
        if (i > 1L) stopifnot(!any(vapply(boxes[seq_len(i - 1L)], intersects, logical(1), a = boxes[[i]])))
        for (j in seq_len(nrow(output$nodes))) {
          node <- output$nodes[j, ]
          stopifnot(!intersects(boxes[[i]], c(node$x - node$width / 2, node$x + node$width / 2,
                                             node$y - node$height / 2, node$y + node$height / 2)))
        }
      }
    }
    grDevices::dev.off()
  }
  cat("All compact path diagram layout checks passed.\n")
}
