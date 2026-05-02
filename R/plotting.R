#' Plot Response Curves for Two Species
#'
#' @param spa Data frame for Species A (x, y)
#' @param spb Data frame for Species B (x, y)
#' @param main Title for the plot
#' @param xlab Label for x axis
#' @param col_a Color for Species A
#' @param col_b Color for Species B
#' @export
plot_response_curves <- function(spa, spb, main = "Response Curves", xlab = "Environmental Gradient", col_a = "black", col_b = "blue") {
  plot(spa$x, spa$y, type = "l", col = col_a, lwd = 2, main = main, xlab = xlab, ylab = "Suitability",
       ylim = c(0, 1), xlim = range(c(spa$x, spb$x)))
  lines(spb$x, spb$y, col = col_b, lwd = 2)

  ov <- beta_overlap(spa, spb)
  if (nrow(ov) > 0) {
    polygon(c(ov$x, rev(ov$x)), c(ov$y, rep(0, nrow(ov))), col = rgb(0.5, 0.5, 0.5, 0.3), border = NA)
  }

  legend("topright", legend = c("Species A", "Species B", "Overlap"), col = c(col_a, col_b, "grey"), lwd = c(2, 2, NA),
         fill = c(NA, NA, rgb(0.5, 0.5, 0.5, 0.3)), border = c(NA, NA, "grey"), bty = "n")
}

#' Plot Niche Divergence Plane
#'
#' @param indices Data frame with exclusivity and dissimilarity columns
#' @param labels Optional vector of labels for the points
#' @param main Title for the plot
#' @export
plot_ndp <- function(indices, labels = NULL, main = "Niche Divergence Plane") {
  plot(indices$exclusivity, indices$dissimilarity, xlim = c(0, 1), ylim = c(0, 1),
       pch = 21, bg = "steelblue", cex = 1.5,
       xlab = "Niche Exclusivity", ylab = "Niche Dissimilarity", main = main)

  abline(h = 0.5, v = 0.5, lty = 2, col = "grey")

  if (!is.null(labels)) {
    text(indices$exclusivity, indices$dissimilarity, labels = labels, pos = 3, cex = 0.8)
  }

  # Add quadrant labels
  text(0.25, 0.25, "Stability", col = "grey", cex = 0.8)
  text(0.75, 0.25, "Expansion", col = "grey", cex = 0.8)
  text(0.25, 0.75, "Shift", col = "grey", cex = 0.8)
  text(0.75, 0.75, "Divergence", col = "grey", cex = 0.8)
}
