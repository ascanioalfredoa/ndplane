#' Determines response overlap between two beta function curves from betaPDF()
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
    ov <-
        spa[spa[, 1] %in% spb[, 1], ] |>
        dplyr::left_join(spb[spb[, 1] %in% spb[, 1], ], by = "x", suffix = c("_a", "_b")) |>
        tidyr::pivot_longer(cols = 2:3, names_to = "source_curve", values_to = "y") |>
        dplyr::group_by(.data$x) |>
        dplyr::summarise(y = min(.data$y))
    ov
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
    1 - (pracma::trapz(ov[[1]], ov[[2]])/pracma::trapz(spa[[1]], spa[[2]]) + pracma::trapz(ov[[1]], ov[[2]])/pracma::trapz(spb[[1]], spb[[2]]))/2
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
    excl <- 1 - (min(c(max(spa[, 1]), max(spb[, 1]))) - max(c(min(spa[, 1]), min(spb[, 1]))))/(max(c(max(spa[, 1]), max(spb[, 1]))) - min(c(min(spa[, 1]), min(spb[, 1]))))
    if(excl > 1) excl <- 1
    excl
}

