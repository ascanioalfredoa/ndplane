#' Thin occurrence records based on pairwise distances
#'
#' Applies spatial thinning to remove points closer than a specified minimum distance. Assumes geographic coordinates (EPSG:4326) unless otherwise specified. Projects to EPSG:4326 if necessary.
#'
#' @param x A matrix, data.frame (lon/lat), or `SpatVector` of coordinates.
#' @param thin.par Minimum pairwise distance in kilometers.
#' @param reps Number of thinning replicates to perform (default = 10).
#'
#' @return A matrix or data.frame of thinned coordinates (longitude, latitude).
#' @export
#' @family thinning
old_thin_records <- function(x, thin.par, reps = 10) {
    if (inherits(x, "SpatVector")) {
        current_crs <- terra::crs(x, describe = TRUE)
        message("Input is a SpatVector. Current CRS: ", current_crs)

        if (terra::crs(x) != "EPSG:4326") {
            message("Projecting input to EPSG:4326 for distance calculations.")
            x <- terra::project(x, "EPSG:4326")
        }

        x <- as.data.frame(terra::geom(x)[, c("x", "y")])
    } else {
        warning("Assuming input coordinates are in EPSG:4326 (WGS84).")
        x <- as.data.frame(x)
        names(x) <- c("x", "y")
    }

    d <- nrow(x)
    results <- vector("list", reps)
    counts <- numeric(reps)

    for (i in 1:reps) {
        sel <- sample(1:d, 1)
        selected <- x[sel, , drop = FALSE]
        to_eval <- x[-sel, , drop = FALSE]

        while (nrow(to_eval) > 0) {
            dist <- terra::distance(terra::vect(to_eval, geom = c("x", "y"), crs = "EPSG:4326"),
                                    terra::vect(selected, geom = c("x", "y"), crs = "EPSG:4326"))
            if (is.matrix(dist)) dist <- apply(dist, 1, min)
            keep <- which(dist >= thin.par * 1000)

            if (length(keep) == 0) break
            sel <- sample(keep, 1)
            selected <- rbind(selected, to_eval[sel, , drop = FALSE])
            to_eval <- to_eval[-sel, , drop = FALSE]
        }

        results[[i]] <- selected
        counts[i] <- nrow(selected)
    }

    return(results[[which.max(counts)]])
}

#' Adaptive grid-based thinning of spatial points
#'
#' Recursively applies spatial thinning by gridding the landscape and thinning dense cells using `thin_records()`.
#' If inputs are not spatial, assumes EPSG:4326 coordinates.
#'
#' @param coords A data.frame or matrix with longitude and latitude columns.
#' @param cellsize Numeric grid cell size in degrees (default = 1).
#' @param min_occ Minimum occurrences per grid cell to apply thinning (default = 10).
#' @param thin_par Minimum distance in kilometers to enforce during thinning.
#' @param reps Number of repetitions for each thinning call (default = 10).
#'
#' @return A data.frame of spatially thinned coordinates.
#' @export
#' @family thinning
old_adaptive_thin <- function(coords, cellsize = 1, min_occ = 10, thin_par = 5, reps = 10) {
    if (!inherits(coords, "sf")) {
        warning("Assuming input coordinates are in EPSG:4326 (WGS84).")
        sf_points <- sf::st_as_sf(coords, coords = c(1, 2), crs = 4326)
    } else {
        current_crs <- sf::st_crs(coords)$input
        message("Input is an sf object. Current CRS: ", current_crs)
        if (sf::st_crs(coords) != sf::st_crs(4326)) {
            message("Projecting input to EPSG:4326 for grid-based processing.")
            coords <- sf::st_transform(coords, 4326)
        }
        sf_points <- coords
    }

    grid <- sf::st_make_grid(sf_points, cellsize = cellsize, what = "polygons")
    intersected <- sf::st_intersects(sf_points, grid)

    keep_all <- list()

    for (i in seq_along(grid)) {
        idx <- which(lengths(intersected) > 0 & vapply(intersected, function(x) i %in% x, logical(1)))
        if (length(idx) >= min_occ && length(idx) < 1000) {
            thinned <- thin_records(sf::st_coordinates(sf_points[idx, ]), thin.par = thin_par, reps = reps)
            keep_all[[length(keep_all) + 1]] <- thinned
        } else if (length(idx) >= 1000) {
            sub_coords <- sf::st_coordinates(sf_points[idx, ])
            recursive <- old_adaptive_thin(sub_coords, cellsize = cellsize / 2, thin_par = thin_par, reps = reps)
            keep_all[[length(keep_all) + 1]] <- recursive
        } else if (length(idx) > 0) {
            keep_all[[length(keep_all) + 1]] <- sf::st_coordinates(sf_points[idx, ])
        }
        colnames(keep_all[[length(keep_all)]]) <- c("x", "y")
    }

    return(do.call(rbind, keep_all))
}


