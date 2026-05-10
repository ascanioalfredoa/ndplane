#' Determines response overlap between two beta function curves
#'
#' @param spa Data frame containing x and y values for Species A response
#' @param spb Data frame containing x and y values for Species B response
#'
#' @return Data frame containing x and y values for the overlap curve between species A and B
#' @export
#'
#' @examples
#' spa <- betaPDF(a = 0, b = 1, alpha = 3, gamma = 1)
#' spb <- betaPDF(a = 0, b = 1, alpha = 1, gamma = 3)
#' beta_overlap(spa, spb)
beta_overlap <- function(spa, spb) {
  # Merge spa and spb by x using base R merge
  merged <- merge(spa, spb, by = "x", suffixes = c("_a", "_b"))

  if (nrow(merged) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0)))
  }

  # Calculate minimum y for each x
  merged$y <- pmin(merged$y_a, merged$y_b)

  # Return data frame with x and y
  merged[, c("x", "y")]
}

#' Calculate Niche Dissimilarity Index for the Niche Divergence Plane
#'
#' @param spa Data frame containing x and y values for Species A response
#' @param spb Data frame containing x and y values for Species B response
#' @param ov Data frame containing x and y values for Overlap between A and B
#'
#' @return Numeric Niche Dissimilarity value
#' @export
#'
#' @examples
#' spa <- betaPDF(a = 0, b = 1, alpha = 3, gamma = 1)
#' spb <- betaPDF(a = 0, b = 1, alpha = 1, gamma = 3)
#' ov <- beta_overlap(spa, spb)
#' niche_diss(spa, spb, ov)
niche_diss <- function(spa, spb, ov) {
  if (nrow(ov) < 2) return(1)

  area_ov <- trapz(ov$x, ov$y)
  area_a <- trapz(spa$x, spa$y)
  area_b <- trapz(spb$x, spb$y)

  1 - (area_ov / area_a + area_ov / area_b) / 2
}

#' Calculate Niche Exclusivity Index for the Niche Divergence Plane
#'
#' @param spa Data frame containing x and y values for Species A response
#' @param spb Data frame containing x and y values for Species B response
#'
#' @return Numeric Niche Exclusivity value
#' @export
#'
#' @examples
#' spa <- betaPDF(a = 0, b = 1, alpha = 3, gamma = 1)
#' spb <- betaPDF(a = 0, b = 1, alpha = 1, gamma = 3)
#' niche_excl(spa, spb)
niche_excl <- function(spa, spb) {
    min_a <- min(spa$x)
    max_a <- max(spa$x)
    min_b <- min(spb$x)
    max_b <- max(spb$x)

    overlap_min <- max(min_a, min_b)
    overlap_max <- min(max_a, max_b)

    total_min <- min(min_a, min_b)
    total_max <- max(max_a, max_b)

    overlap_width <- max(0, overlap_max - overlap_min)
    total_width <- total_max - total_min

    excl <- 1 - (overlap_width / total_width)
    if(excl > 1) excl <- 1
    if(excl < 0) excl <- 0
    excl
}
