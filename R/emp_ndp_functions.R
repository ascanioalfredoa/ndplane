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

#' Prepare background points with Target Group Species selection and thinning
#'
#' @param all_occurrences list List of all salamander occurrence spatial objects
#' @param thinning_params list Parameters for thinning background points (optional)
#'
#' @return list containing processed background points
#' @export
prepare_background_points <- function(all_occurrences, thinning_params=NULL) {
  # TODO: combine salamander occurrences to generate TGS
  # TODO: thin background points as needed
  # return processed background points
  list(background_points = NULL)
}

#' Load environmental rasters, extract raster values for occurrences and background, remove incomplete cases
#'
#' @param env_raster_paths character Vector of paths to environmental raster files
#' @param occurrence_points SpatialPoints or SpatialPointsDataFrame Occurrence points
#' @param background_points SpatialPoints or SpatialPointsDataFrame Background points
#'
#' @return list with occurrence_data and background_data data frames
#' @export
load_and_extract_environmental_data <- function(env_raster_paths, occurrence_points, background_points) {
  # TODO: load rasters
  # TODO: extract raster values for points
  # TODO: filter out points with missing data
  list(occurrence_data = NULL, background_data = NULL)
}

#' Run maxent analysis for one environmental variable
#'
#' @param env_variable character Name of environmental variable to use
#' @param species_occurrences list List of species-specific occurrence data frames
#' @param species_backgrounds list List of species-specific background data frames
#'
#' @return list containing ind_res and plot_res result objects
#' @export
run_maxent_analysis_for_variable <- function(env_variable, species_occurrences, species_backgrounds) {
  # TODO: create ENMTools species objects
  # TODO: create maxent models
  # TODO: calculate Schoener's D, Hellinger's I, and rank correlation
  # TODO: run maxent_ndp function
  # separate and return analysis results
  list(ind_res = NULL, plot_res = NULL)
}

#' Save results (indices and plotting data) to CSV or RDS
#'
#' @param ind_res data.frame Results containing indices
#' @param plot_res data.frame Results containing plotting information
#' @param ind_res_filepath character File path to save indices results
#' @param plot_res_filepath character File path to save plotting results
#'
#' @return NULL
#' @export
save_results <- function(ind_res, plot_res, ind_res_filepath, plot_res_filepath) {
  # TODO: write ind_res and plot_res to files
}
