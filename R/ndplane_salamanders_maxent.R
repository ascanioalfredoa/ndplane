#' Load Vector Data
#'
#' Loads spatial vector data from a shapefile.
#' @param filepath Character. File path to the shapefile.
#' @return Spatial vector object.
#' @import terra
#' @export
load_vector_data <- function(filepath) {
  terra::vect(filepath)
}

#' Create Study Area Buffers
#'
#' Creates buffered study areas around occurrences and intersects with a basemap.
#' @param occurrences Spatial vector of occurrence points.
#' @param basemap Spatial vector background map.
#' @param buffer_distance Numeric. Buffer distance in map units.
#' @return Buffered and intersected spatial vector study area.
#' @import terra
#' @export
create_study_area <- function(occurrences, basemap, buffer_distance) {
  sa <- terra::buffer(occurrences, buffer_distance)
  sa <- terra::aggregate(sa)
  sa <- terra::intersect(sa, basemap)
  sa
}

#' Thin Records
#'
#' Wrapper to thinning function (expects thin_records to be available in namespace).
#' @param records Spatial vector of records to thin.
#' #@param thin_par Numeric. Parameter for thinning distance.
#' @param reps Integer. Number of repetitions.
#' @return Thinned spatial vector.
#' @export
thin_occurrence_records <- function(records, thin_par = 10, reps = 10) {
  thin_records(records, thin.par = thin_par, reps = reps)
}

#' Prepare Environmental Data
#'
#' Load and crop environmental raster layers to study area.
#' @param raster_paths Character vector of file paths to raster layers.
#' @param study_area Spatial vector of study area.
#' @return SpatRaster of cropped environmental data.
#' @import terra
#' @export
prepare_env_data <- function(raster_paths, study_area) {
  env_stack <- terra::rast(raster_paths)
  env_crop <- terra::crop(env_stack, study_area)
  env_mask <- terra::mask(env_crop, study_area)
  env_mask
}

#' Extract Environmental Values
#'
#' Extract environmental values at spatial points.
#' @param env_data SpatRaster of environmental data.
#' @param points Spatial vector of points.
#' @return Data frame of extracted environmental values.
#' @import terra
#' @export
extract_env_values <- function(env_data, points) {
  terra::extract(env_data, points)
}

#' Setup Presence-Absence Data
#'
#' Combine presence and background points for Maxent input.
#' @param presence Spatial vector of presence points.
#' #@param background Spatial vector of background points.
#' @param env_data Environmental data extracted at points.
#' @return Data frame of presence-absence and environmental variables.
#' @export
setup_pa_data <- function(presence, background, env_data) {
  pa <- c(rep(1, length(presence)), rep(0, length(background)))
  data <- rbind(presence, background)
  data_extracted <- terra::extract(env_data, data)[, -1, drop = FALSE]
  combined <- data.frame(pa = pa, data_extracted)
  combined <- combined[complete.cases(combined),]
  combined
}

#' Fit Univariate Maxent Model
#'
#' Fit maxnet model for a single environmental variable.
#' @param pa_data Data frame with presence-absence and variables.
#' @param env_var Character name of environmental variable.
#' @return maxnet model object.
#' @import maxnet
#' @export
fit_maxent_univariate <- function(pa_data, env_var) {
  maxnet::maxnet(pa_data$pa, pa_data[[env_var]], maxnet::maxnet.formula( ~ ., pa_data))
}

#' Calculate Niche Divergence Indices
#'
#' Calculate niche similarity/dissimilarity metrics from maxent models.
#' @param model1 maxnet model for species 1.
#' @param model2 maxnet model for species 2.
#' #@param env RasterLayer or numeric vector of environmental data
#' @return List of niche indices and response curves.
#' @export
calculate_maxent_ndp <- function(model1, model2) {
  # Placeholder: Implement using ENMTools or custom niche metrics
  list(indices = NULL, response_curves = NULL)
}

#' Run Full Maxent NDP Analysis
#'
#' Runs the full pipeline: load data, create study areas, thin data, run maxent per variable, calculate niche indices.
#' @param presence_species1 Path to presence shapefile for species 1.
#' @param presence_species2 Path to presence shapefile for species 2.
#' @param background_species1 Path to background points shapefile for species 1.
#' @param background_species2 Path to background points shapefile for species 2.
#' @param basemap_path Path to basemap shapefile.
#' @param raster_paths Character vector of environmental raster file paths.
#' @param buffer_distance Numeric buffer size for study areas.
#' @return List with niche indices, response curves, and other data.
#' @export
run_maxent_ndp_pipeline <- function(presence_species1, presence_species2,
                                    background_species1, background_species2,
                                    basemap_path, raster_paths, buffer_distance = 1e5) {
  
  # Load spatial data
  sp1 <- load_vector_data(presence_species1)
  sp2 <- load_vector_data(presence_species2)
  bg1 <- load_vector_data(background_species1)
  bg2 <- load_vector_data(background_species2)
  basemap <- load_vector_data(basemap_path)
  
  # Create study areas
  sa1 <- create_study_area(sp1, basemap, buffer_distance)
  sa2 <- create_study_area(sp2, basemap, buffer_distance)
  all_sa <- terra::union(sa1, sa2) %>% terra::aggregate()
  
  # Prepare environmental data
  env_data <- prepare_env_data(raster_paths, all_sa)
  
  # Thin background points (if needed)
  bg1_thin <- thin_occurrence_records(bg1)
  bg2_thin <- thin_occurrence_records(bg2)
  
  # Loop through raster layers
  niche_indices_list <- list()
  response_curves_list <- list()
  for(i in seq_along(names(env_data))) {
    layer_name <- names(env_data)[i]
    
    # Setup presence-absence data
    pa_data1 <- setup_pa_data(sp1, bg1_thin, env_data[[i]])
    pa_data2 <- setup_pa_data(sp2, bg2_thin, env_data[[i]])
    
    # Fit maxent models
    model1 <- fit_maxent_univariate(pa_data1, layer_name)
    model2 <- fit_maxent_univariate(pa_data2, layer_name)
    
    # Calculate niche divergence indices
    ndp_res <- calculate_maxent_ndp(model1, model2)
    
    niche_indices_list[[layer_name]] <- ndp_res$indices
    response_curves_list[[layer_name]] <- ndp_res$response_curves
  }
  
  list(niche_indices = niche_indices_list, response_curves = response_curves_list)
}
