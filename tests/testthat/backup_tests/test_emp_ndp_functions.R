# Create example SpatVectors for testing
example_crs <- "EPSG:4326"
point1 <- terra::vect(matrix(c(-100, 40), ncol=2), type="points", crs=example_crs)
point2 <- terra::vect(matrix(c(-101, 41), ncol=2), type="points", crs=example_crs)

# Create another SpatVector with different CRS
example_crs2 <- "EPSG:3857"
point3 <- terra::project(point1, example_crs2)

# Write a test for same_crs
test_that("same_crs correctly detects matching CRS", {
  expect_true(same_crs(list(point1, point2)))
  expect_false(same_crs(list(point1, point3)))
})

# Write a test for read_point_shapefiles
# For this, we create temporary shapefiles

shp_temp1 <- tempfile(fileext = ".shp")
shp_temp2 <- tempfile(fileext = ".shp")
terra::writeVector(point1, shp_temp1, overwrite=TRUE)
terra::writeVector(point2, shp_temp2, overwrite=TRUE)



# The test for read_point_shapefiles should read these and confirm they are SpatVectors with points


test_that("read_point_shapefiles reads shapefiles and validates geometry", {
  svs <- read_point_shapefiles(c(shp_temp1, shp_temp2))
  expect_type(svs, "list")
  expect_true(all(sapply(svs, inherits, what = "SpatVector")))
  # All should be points
  expect_true(all(sapply(svs, function(sv) terra::geomtype(sv) %in% c("points", "point"))))
})

# Create example data for testing create_multi_sa
# Create two sets of points representing two species
species1 <- terra::vect(data.frame(x = c(1,2,3), y = c(1,2,3)), geom = c("x", "y"), crs = "EPSG:4326")
species2 <- terra::vect(data.frame(x = c(4,5,6), y = c(4,5,6)), geom = c("x", "y"), crs = "EPSG:4326")

occ_list <- list(species1 = species1, species2 = species2)

test_that("create_multi_sa works on lists of spatvectors", {
  # Use create_multi_sa to create study areas
  sa_list <- create_multi_sa(occ_list, buffer = 1000)
  
  # Some expectations
  expect_true(is.list(sa_list))
  expect_equal(length(sa_list), 2)
  expect_equal(names(sa_list), c("species1", "species2"))
  expect_true(all(sapply(sa_list, function(x) inherits(x, "SpatVector"))))
})  

# Optionally, test create_multi_sa with a single SpatVector with taxa_field
species1$taxa <- "species1"
species2$taxa <- "species2"
merged_sv <- rbind(species1, species2)

test_that("create_multi_sa works on a single object with a field representing taxa", {
  sa_split <- create_multi_sa(merged_sv, buffer = 1000, taxa_field = "taxa")
  expect_true(is.list(sa_split))
  expect_equal(length(sa_split), 2)
  expect_true(all(sapply(sa_split, function(x) inherits(x, "SpatVector"))))
})

# Test create_multi_sa with single SpatVector no taxa_field (should call create_sa once)
test_that("create_multi_sa with single SpatVector calls create_sa", {
  result <- create_multi_sa(point1, buffer = 100, bg_area = NULL, taxa_field = NULL)
  expect_true(inherits(result, "SpatVector")) 
})

# Test create_multi_sa with list of SpatVectors
sp_list <- list(point1 = point1, point2 = point2)
test_that("create_multi_sa with list of SpatVectors calls create_sa per species", {
  result <- create_multi_sa(sp_list, buffer=100, bg_area = NULL)
  expect_type(result, "list")
  expect_length(result, 2)
})
