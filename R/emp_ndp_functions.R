#' Create study area proxy using buffers around occurrences
#'
#' @param occs list List of species occurrence spatial objects
#' @param buffer numeric Buffer distance in meters (default 100000 = 100 km)
#' @param bg_area numeric Buffer distance in meters (default 100000 = 100 km)
#'
#' @return list containing study area proxy based on buffer around occurrence points
#' @export
create_sa <- function(occs, buffer = 1e5, bg_area = NULL) {
  # Validate occs input
  if(!inherits(occs, "SpatVector")) stop("'occs' must be a terra SpatVector object.")
  
  # Create buffer around occurrences using buffer parameter
  sa <- terra::buffer(occs, buffer)
  # Aggregate and merge overlapping buffers
  sa <- terra::aggregate(sa)
  
  if(!is.null(bg_area)) {
    # If bg_area is a path to shapefile, read it
    if(is.character(bg_area)) {
      bg_area <- tryCatch(terra::vect(bg_area), 
                          error = function(e) stop("Could not read 'bg_area' from provided path."))
    }
    # Check if bg_area is polygon SpatVector
    if(!terra::is.polygons(bg_area)) {
      stop("'bg_area' should be a polygon SpatVector or a path to polygon shapefile.")
    }
    # Intersect buffers with background polygon
    sa <- terra::intersect(sa, bg_area)
    if(nrow(sa) == 0) warning("Intersection result of buffer and bg_area is empty.")
  }
  
  return(sa)
}

#' Check if all SpatVectors have the same CRS
#'
#' @param svs list of SpatVectors
#' @return logical TRUE if all have the same CRS, FALSE otherwise
  same_crs <- function(svs) {
    crs_list <- sapply(svs, function(x) terra::crs(x))
    all(crs_list == crs_list[1])
  }
    
#' Read and validate list of shapefiles containing points
#'
#' @param file_paths character vector of file paths to shapefiles
#' @return list of SpatVectors with point geometries
read_point_shapefiles <- function(file_paths) {
  taxa_occs <- lapply(file_paths, function(fp) {
    sv <- terra::vect(fp)
    if(!terra::geomtype(sv) %in% c("points", "point")) {
      stop(paste0("File '", fp, "' does not contain point geometries."))
    }
    sv
  })
  taxa_occs
}

#' Create study areas for multiple species occurrences
#'
#' @param occs_list list or SpatVector List of species occurrence spatial objects or a single SpatVector with a taxa/species identifier column
#' @param buffer numeric Buffer distance in meters (default 100000 = 100 km)
#' @param bg_area numeric Buffer distance in meters (default 100000 = 100 km)
#' @param taxa_field character Name of the taxa/species identifier field in single SpatVector. Default NULL.
#' @param desired_crs character EPSG code or CRS string to reproject all SpatVectors to (optional)
#'
#' @return named list of study area proxies for each species
#' @export
create_multi_sa <- function(occs_list, buffer = 1e5, bg_area = NULL, taxa_field = NULL, desired_crs = NULL) {

  # If occs_list is a single SpatVector
  if(inherits(occs_list, "SpatVector")) {
    # If taxa_field is provided, split by taxa
    if(!is.null(taxa_field)) {
      if(!taxa_field %in% names(occs_list)) {
        stop(paste0("Field '", taxa_field, "' not found in the SpatVector attributes."))
      }
      taxa_values <- unique(occs_list[[taxa_field]]) |> unlist()
      # Create named list of SpatVectors by taxa
      taxa_occs <- lapply(taxa_values, function(tv) {
        occs_list[occs_list[[taxa_field]] == tv, ]
      })
      names(taxa_occs) <- taxa_values
    } else {
      # No taxa_field provided, run regular create_sa for the single SpatVector
      sa <- create_sa(occs_list, buffer = buffer, bg_area = bg_area)
      return(sa)
    }
  } else if(is.list(occs_list)) {
    # Check if list elements are all SpatVectors
    if(all(sapply(occs_list, inherits, what = "SpatVector"))) {
      taxa_occs <- occs_list
    } else if(all(sapply(occs_list, function(x) is.character(x) && file.exists(x)))) {
      # Assume list of file paths, read the shapefiles
      taxa_occs <- read_point_shapefiles(occs_list)
    } else {
      stop("Input list must consist of SpatVectors or paths to shapefiles containing points.")
    }

    # If list has no names, assign default names
    if(is.null(names(taxa_occs))) {
      names(taxa_occs) <- paste0("species_", seq_along(taxa_occs))
    }

    # Check projections, reproject as needed
    if(!same_crs(taxa_occs)) {
      if(is.null(desired_crs)) {
        # Project all to CRS of first SpatVector
        target_crs <- terra::crs(taxa_occs[[1]])
      } else {
        # Use user-provided CRS
        target_crs <- desired_crs
      }
      taxa_occs <- lapply(taxa_occs, function(sv) {
        if(terra::crs(sv) != target_crs) {
          terra::project(sv, target_crs)
        } else {
          sv
        }
      })
    }
  } else {
    # If not a SpatVector or list, run regular create_sa
    sa <- create_sa(occs_list, buffer = buffer, bg_area = bg_area)
    return(sa)
  }

  # Iterate over taxa_occs and create study areas
  sa_list <- lapply(taxa_occs, function(sp_occs) {
    create_sa(sp_occs, buffer = buffer, bg_area = bg_area)
  })

  return(sa_list)
}

#' Load environmental rasters, extract raster values for occurrences and background, remove incomplete cases
#'
#' @param ras character Vector of paths to environmental raster files or SpatRaster
#' @param occ SpatVector or object convertible to SpatVector Occurrence points
#' @param bgp SpatVector or object convertible to SpatVector Background points
#'
#' @return list with occurrence_data and background_data data frames
#' @export
extract_environmental_data <- function(ras, occ, bgp) {

  # Load or use provided SpatRaster
  if (inherits(ras, "SpatRaster")) {
    env_rasters <- ras
  } else {
    env_rasters <- terra::rast(ras)
  }

  # Load or use provided SpatVector for occurrences
  if (inherits(occ, "SpatVector")) {
    occ <- occ
  } else {
    occ <- terra::vect(occ)
  }

  # Load or use provided SpatVector for background points
  if (inherits(bgp, "SpatVector")) {
    bgp <- bgp
  } else {
    bgp <- terra::vect(bgp)
  }

  # Create vector of presence (1) and absence (0)
  pa <- c(rep(1, length(occ)), rep(0, length(bgp)))

  # Combine presence absence datasets
  combined_points <- rbind(occ, bgp)

  # Extract environmental values using the new set of points
  envdat <- terra::extract(env_rasters, combined_points)[, -1, drop = FALSE]

  # Add pa vector as a column to extracted values
  data <- cbind(pa = pa, envdat)

  # Clean NAs
  data <- data[complete.cases(data), ]

  # Separate pa and envdat again
  pa <- data$pa
  envdat <- data[, -1, drop = FALSE]

  list(occurrence_data = list(pa = pa[pa == 1], env = envdat[pa == 1, , drop = FALSE]),
       background_data = list(pa = pa[pa == 0], env = envdat[pa == 0, , drop = FALSE]))
}

#' Calculate niche divergence plane indices based on maxent response curves
#'
#' @param mx_sp1 maxent model object 1 (EMNTools::enmtools.maxent output)
#' @param mx_sp2 maxent model object 2 (EMNTools::enmtools.maxent output)
#' @param cut_off numeric Proportion cutoff to trim response curve tails (default 0.95)
#'
#' @return list containing:
#'   - indices: data.frame with dissimilarity and exclusivity indices
#'   - response_curves: combined trimmed response curves of the two species
#'
#' @importFrom pracma trapz
#' @importFrom zoo na.approx
#' @export
maxent_ndp <- function(mx_sp1, mx_sp2, cut_off = 0.95) {
  # Save variable name
  var <- unique(names(mx_sp1$response.plots)[[1]], names(mx_sp1$response.plots)[[1]])

  # Save species responses in different objects
  A <- mx_sp1$response.plots[[1]]$data |> dplyr::filter(source == "Suitability")
  B <- mx_sp2$response.plots[[1]]$data |> dplyr::filter(source == "Suitability")

  A$source <- mx_sp1$species.name
  B$source <- mx_sp2$species.name

  # Recalculate original probability value to assess quantiles
  A$cs <- cumsum(A$value / sum(A$value))
  B$cs <- cumsum(B$value / sum(B$value))

  # Cut responses given the cut_off threshold
  lco <- (1 - cut_off) / 2
  uco <- 1 - ((1 - cut_off) / 2)

  A <- A[A$cs >= lco & A$cs <= uco, ]
  B <- B[B$cs >= lco & B$cs <= uco, ]

  #### Integration by determining the lower function ####
  # Where the values of X are shared between distributions
  # And where the values of Y are the lowest between the two possibilities

  C <-
    A[A[, 1] >= min(B[, 1]) & A[, 1] <= max(B[, 1]), 1:3] |>
    dplyr::full_join(B[B[, 1] >= min(A[, 1]) & B[, 1] <= max(A[, 1]), 1:3], by = "layer") |>
    dplyr::arrange(layer) |>
    dplyr::mutate(value.x = zoo::na.approx(value.x, x = layer, na.rm = FALSE)) |>
    dplyr::mutate(value.y = zoo::na.approx(value.y, x = layer, na.rm = FALSE)) |>
    dplyr::select(1, 2, 4) |>
    dplyr::filter(!is.na(value.x) & !is.na(value.y)) |>
    tidyr::pivot_longer(cols = 2:3, names_to = "source_curve", values_to = "y") |>
    dplyr::group_by(layer) |>
    dplyr::summarise(y = min(y))

  # Calculate area
  NicheDiss <- 1 - (pracma::trapz(C[[1]], C[[2]]) / pracma::trapz(A[[1]], A[[2]]) +
                     pracma::trapz(C[[1]], C[[2]]) / pracma::trapz(B[[1]], B[[2]])) / 2

  NicheEx <- 1 - (min(c(max(A[[1]]), max(B[[1]]))) - max(c(min(A[[1]]), min(B[[1]])))) /
                   (max(c(max(A[[1]]), max(B[[1]]))) - min(c(min(A[[1]]), min(B[[1]]))))
  if(NicheEx > 1) NicheEx <- 1

  res <- list()
  res[["indices"]] <- data.frame(dissimilarity = NicheDiss, exclusivity = NicheEx, var = var)
  res[["response_curves"]] <- rbind(A, B)
  res[["response_curves"]]$var <- var

  return(res)
}

#' Run maxent analysis for one or more environmental variables
#'
#' @param envdata SpatRaster Environmental raster data
#' @param species_occurrences list List of species-specific occurrence SpatVectors
#' @param species_backgrounds list List of species-specific background SpatVectors
#' @param study_areas list List of study areas (SpatVectors) corresponding to species
#' @param variables character or numeric "all" (default) or a subset of variable names or indices
#'
#' @return list containing ind_res and plot_res result objects
#' @export
run_maxent_ndp <- function(
  envdata,
  species_occurrences,
  species_backgrounds,
  study_areas,
  variables = "all"
) {

  # Check if variables is all or subset
  if (is.character(variables) && length(variables) > 1) {
    # vector of names (character)
    if (!all(variables %in% names(envdata))) {
      stop("Some variable names are not in envdata")
    }
    var_indices <- which(names(envdata) %in% variables)
  } else if (is.character(variables) && length(variables) == 1) {
    if (variables == "all") {
      var_indices <- 1:terra::nlyr(envdata)  # number of layers in envdata raster
    } else if (variables %in% names(envdata)) {
      var_indices <- which(names(envdata) == variables)
    } else {
      stop("Variable '", variables, "' not found in envdata")
    }
  } else if (is.numeric(variables)) {
    if (all(variables >= 1 & variables <= terra::nlyr(envdata))) {
      var_indices <- variables
    } else {
      stop("Numeric variable indices out of range")
    }
  } else {
    stop("Invalid type for variables parameter")
  }

  ind_res <- NULL
  plot_res <- NULL

  for (i in var_indices) {
    # create species objects
    sp1_et <- ENMTools::enmtools.species(
      range = terra::mask(envdata[[i]], study_areas[[1]]),
      presence.points = species_occurrences[[1]],
      background.points = species_backgrounds[[1]],
      species.name = names(species_occurrences)[1]
    )

    sp2_et <- ENMTools::enmtools.species(
      range = terra::mask(envdata[[i]], study_areas[[2]]),
      presence.points = species_occurrences[[2]],
      background.points = species_backgrounds[[2]],
      species.name = names(species_occurrences)[2]
    )
    # run maxent models
    sp1_maxent <- ENMTools::enmtools.maxent(sp1_et, env = envdata[[i]], bg.source = "points")
    sp2_maxent <- ENMTools::enmtools.maxent(sp2_et, env = envdata[[i]], bg.source = "points")
    # calculate Schoener's D and Hellinger's I
    mx_sim <- ENMTools::raster.overlap(sp1_maxent, sp2_maxent) |>
      (\(x) do.call(rbind, x))()

    # calculate niche divergence plane indices
    ndp_sim <- maxent_ndp(mx_sp1 = sp1_maxent, mx_sp2 = sp2_maxent, cut_off = 0.95)

    ind_res <- rbind(ind_res, cbind(ndp_sim$indices, t(mx_sim)))
    plot_res <- rbind(plot_res, ndp_sim$response_curves)
}

  list(ind_res = ind_res, plot_res = plot_res)
}
