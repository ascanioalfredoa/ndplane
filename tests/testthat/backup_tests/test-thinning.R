test_that("thin_records returns expected structure and CRS assumption warning", {
    skip_if_not_installed("terra")

    # Simulate simple points
    set.seed(42)
    pts <- data.frame(
        lon = runif(20, -85, -80),
        lat = runif(20, 35, 40)
    )

    expect_warning({
        thinned <- thin_records(pts, thin.par = 10, reps = 5)
    }, regexp = "Assuming input coordinates are in EPSG:4326")

    expect_type(thinned, "list")
    expect_equal(ncol(thinned), 2)
    expect_true(nrow(thinned) <= 20)
})

test_that("thin_records handles SpatVector input with CRS message", {
    skip_if_not_installed("terra")

    pts <- data.frame(
        lon = runif(15, -90, -85),
        lat = runif(15, 30, 35)
    )
    spv <- terra::vect(pts, crs = "EPSG:4326")

    expect_message({
        res <- thin_records(spv, thin.par = 5, reps = 3)
    }, regexp = "Input is a SpatVector")
    expect_type(res, "list")
    expect_equal(ncol(res), 2)
})

test_that("adaptive_thin runs and returns expected output", {
    skip_if_not_installed("sf")

    pts <- data.frame(
        lon = runif(50, -120, -115),
        lat = runif(50, 45, 48)
    )

    expect_warning({
        res <- adaptive_thin(pts, cellsize = 1, min_occ = 5, thin_par = 10, reps = 3)
    }, regexp = "Assuming input coordinates are in EPSG:4326")

    expect_type(res, "list")
    expect_equal(ncol(res), 2)
    expect_true(nrow(res) <= 50)
})
