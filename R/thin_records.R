#' Thin spatial records by minimum distance
#'
#' Performs spatial thinning on point data to ensure a minimum separation between
#' records. The algorithm iteratively removes the record with the highest number
#' of neighbors within the threshold, repeats for a number of replicates, and
#' returns the best solution.
#'
#' @param data A data.frame, tibble, matrix, or terra::SpatVector containing
#'   point coordinates. For tabular inputs, columns defined by `lonlat` must be
#'   present.
#' @param lonlat Character vector of length 2 giving the longitude and latitude
#'   column names for tabular inputs. Ignored when `data` is a SpatVector.
#' @param thin.par Numeric. Minimum allowed distance between any two records. The
#'   interpretation of units follows terra::distance (great-circle distance for
#'   lon/lat coordinates; kilometers when `unit = "km"` is used).
#' @param reps Integer. Number of replicate thinnings to attempt; the best
#'   solution (maximizing retained records, then mean pairwise distance) is
#'   returned.
#'
#' @return An object of the same class as `data` containing only the retained
#'   records after thinning.
#'
#' @details Distances are computed with terra::distance. For data.frames/tibbles
#'   the coordinates are assumed to be longitude/latitude in decimal degrees
#'   (great-circle distance). For matrices and SpatVectors distances are computed
#'   in kilometers in this implementation.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' pts <- data.frame(lon = runif(50, -80, -79), lat = runif(50, 8, 9))
#' thn <- thin_records(pts, thin.par = 10, reps = 5)
#' nrow(pts); nrow(thn)
#' }
#' @importFrom dplyr select
#' @importFrom terra distance
#' @export
thin_records <- function(data, lonlat = c("lon", "lat"), thin.par, reps = 10) {
    
    if(any(class(data) == "data.frame") | any(class(data) == "tibble")) {
        spatdat <- data |> dplyr::select(lonlat[1], lonlat[2]) |> as.matrix()
        distmat <- terra::distance(spatdat, lonlat = T) |> 
            as.matrix()
        
    } else if(class(data) == "matrix") {
        spatdat <- data |> dplyr::select(lonlat[1], lonlat[2])
        distmat <- terra::distance(spatdat, lonlat = T, unit = "km") |> 
            as.matrix()
        
    } else if(class(data) == "SpatVector"){
        spatdat <- data
        distmat <- terra::distance(spatdat, unit = "km") |> 
            as.matrix()
    }
    
    res <- 0
    
    for(i in 1:reps) {
        distlogi <- distmat < thin.par
        diag(distlogi) <- F
        sumvec <- rowSums(distlogi)
        keepvec <- sumvec >= 0
        
        
        while (any(distlogi) && sum(keepvec) > 1) {
            RemoveRec <- which(sumvec == max(sumvec))
            
            if (length(RemoveRec) > 1 && any(distlogi[RemoveRec, RemoveRec])) {
                RemoveRec <- sample(RemoveRec, 1)
            }
            
            distlogi[RemoveRec, ] <- FALSE
            distlogi[, RemoveRec] <- FALSE
            keepvec[RemoveRec] <- FALSE
            
            sumvec <- rowSums(distlogi)
            
        }
        
        if(sum(keepvec) > sum(res)) {
            res <- keepvec
        } else if(i > 0 && sum(keepvec) > sum(res)) {
            res <- list(res, keepvec)[[
                which.max(c(mean(distmat[res, res]),
                            mean(distmat[res, res])))
            ]]
        }
    }
    data[res, ]
}

library(sf)

#> Linking to GEOS 3.6.2, GDAL 2.2.3, proj.4 4.9.3

#' Intersect points with a regular grid and index grid cells
#'
#' Builds a regular grid over the extent of the input points, intersects points
#' with the grid, and assigns a sequential cell index.
#'
#' @param p A terra::SpatVector of points.
#' @param cellsize Numeric. Grid cell size (units follow the CRS of `p`).
#'
#' @return A terra::SpatVector resulting from the intersection, with an `Index`
#'   attribute identifying the grid cell for each point.
#'
#' @examples
#' \dontrun{
#' # p is a SpatVector of points
#' A <- gridify_int(p, cellsize = 10)
#' }
#'
#' @importFrom sf st_as_sf st_make_grid
#' @importFrom terra vect intersect
#' @export
gridify_int <- function(p, cellsize = 10) {
    sfp <- sf::st_as_sf(p[, 0])
    grd <- sf::st_make_grid(sfp, cellsize = cellsize)
    grd_v <- terra::vect(grd)
    grd_v$Index <- 1:length(grd_v)
    terra::intersect(p, grd_v)
}

#' Create a grid-based index ensuring a maximum occupancy per cell
#'
#' Applies gridify_int with a starting cell size (computed from the extent if
#' not supplied) and recursively halves the cell size until all grid cells have
#' 1000 or fewer points.
#'
#' @param p A terra::SpatVector of points.
#' @param cs Numeric. Initial grid cell size. If NULL, a default based on
#'   `terra::ext(p)` is used.
#'
#' @return A terra::SpatVector with an `Index` attribute assigning each point to
#'   a grid cell.
#'
#' @examples
#' \dontrun{
#' A <- gridify(p)
#' }
#'
#' @importFrom terra ext
#' @export
gridify <- function(p, cs = NULL) {
    if(is.null(cs)) {
        p_ext <- terra::ext(p)
        cs <- round(min(abs(p_ext[1] - p_ext[2]), abs(p_ext[3] - p_ext[4])), 3)/2
    }
    A <- gridify_int(p, cellsize = cs)
    print(sort(table(A$Index), decreasing = T))
    while(any(table(A$Index) > 1000)) {
        cs <- cs/2
        A <- gridify_int(p, cellsize = cs)
        print(cs)
        print(sort(table(A$Index), decreasing = T))
    }
    A
}

#' Grid-based recursive thinning of spatial points
#'
#' Combines grid partitioning and distance-based thinning. Points in cells with
#' up to 1000 records are thinned using thin_records; cells with more than 1000
#' records are further subdivided until the threshold is met, then thinned. The
#' result is thinned again to enforce the minimum distance globally.
#'
#' @param p A terra::SpatVector of points.
#' @param cs Numeric. Grid cell size used for subdivision. If NULL, a default is
#'   computed from `terra::ext(p)`.
#' @param thin.par Numeric. Minimum allowed distance between any two retained
#'   records (see thin_records).
#' @param reps Integer. Number of replicate thinnings to attempt within each
#'   cell.
#'
#' @return A terra::SpatVector of points after recursive grid partitioning and
#'   distance-based thinning.
#'
#' @examples
#' \dontrun{
#' B <- gridify_thin(p, thin.par = 10)
#' }
#'
#' @importFrom terra ext
#' @export
gridify_thin <- function(p, cs = NULL, thin.par = 10, reps = 10) {
    if(is.null(cs)) {
        p_ext <- ext(p)
        cs <- round(min(abs(p_ext[1] - p_ext[2]), abs(p_ext[3] - p_ext[4])), 3)/2
    }
    A <- gridify_int(p, cellsize = cs)
    print("gridify_int done")
    print(table(A$Index))
    # Split gridcells into two - less equal or greater than 1000 #
    # Gridcells with less-equal than 1000 go through thinning
    # Gridcells with more than 1000 go through more gridify
    # Points at every step get added to a resulting vector and removed from p
    thin_index <- names(which(table(A$Index) <= 1000))
    print("thin_index done")
    res <- p[0, ]
    for(i in thin_index) {
        res <- rbind(res, thin_records(A[A$Index == i, 0], thin.par = thin.par, reps = reps))
        print(i)
    }
    print("thinning done")
    
    grid_index <- names(which(!table(A$Index) <= 1000))
    for(i in grid_index) {
        p2 <- A[A$Index == i, 0]
        res <- rbind(res, gridify_thin(p2, cs = cs/2, thin.par = thin.par, reps = reps))
    }
    
    res <- thin_records(res, thin.par = thin.par)
}

#B <- gridify_thin(p, thin.par = 10)
#B <- gridify(p)
#D <- thin_records(B, thin.par = 10)
#list(unique(B$Index))
#thin_records(B[B$Index %in% unique(B$Index)[1], 0], thin.par = 10)


#D <- lapply(unique(B$Index), function(x) thin_records(B[B$Index %in% unique(B$Index)[x], 0], thin.par = 10))
#points(occs)