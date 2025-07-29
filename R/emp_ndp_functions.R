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
    if(terra::ngeom(sa) == 0) warning("Intersection result of buffer and bg_area is empty.")
  }
  
  return(sa)
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
